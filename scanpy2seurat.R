library(reticulate)
library(anndata)
library(Seurat)
library(Matrix)
library(optparse)

#!/usr/bin/env Rscript

py_require(c("anndata"))

# ---- Functions ----
dgR_flip <- function(dgr){
  #' Helper function that converts scanpy/anndata dgRs to Seurat dgCs
  #' 
  #' Seurat stores cell information in cell x gene format as dgCs
  #' AnnData stores cell information in gene x cell format as dgRs
  #' This function helps convert between them
  #' 
  #' @param dgr input dgR sparse matrix
  #' 
  return(as(as(Matrix::t(dgr), "CsparseMatrix"), "dgCMatrix"))
}


GetScanpyMtx <- function(adata,
                         layer = NULL){
  #' Read information from anndata and flip
  #' 
  #' @param adata annData object, or path to it
  #' @param layer layer containing counts. Default: Get X.
  
  if (is.character(adata)) {
    adata <- read_h5ad(adata)
  }
  
  # Get counts and transform to match SeuratObject conventions
  if( is.null(layer) || (layer == "X") ){
    mtx <- dgR_flip(adata$X)
  } else if (startsWith(layer, "raw")){
    mtx <- dgR_flip(adata$raw$X)
  } else {
    mtx <- dgR_flip(adata$layers[layer])
  }
  return(mtx)
}

Scanpy2Seurat <- function(h5ad,
                          counts = "X",
                          data_counts = "raw",
                          idents = "sample",
                          min_obs= 500,
                          min_var= 300,
                          assay = "RNA",
                          project = "Scanpy",
                          export = TRUE){
#' Convert Scanpy H5ADs to SeuratObject.
#' Retain as much information (ex: HVGs, embeddings) as possible.
#' 
#' @param h5ad path to h5ad file, or h5ad_file
#' @param counts layer containing raw counts
#' @param idents column in adata[["obs"]] to use as main clustering factor
#' @param min_obs Include features detected in at least this many cells.
#' @param min_var Include cells with at least this many features.
#' @param assay Assay name for the final SeuratObject
#' @param project Project name for the final SeuratObject
#' @param export Save as RDS
  
  # Unpack adata information
  if (is.character(h5ad)) {
    adata <- read_h5ad(h5ad)
  }
  
  mtx <- GetScanpyMatrix(adata, counts)
  X <- GetScanpyMatrix(adata, data_counts)
  
  # Save obs and var matrices
  s_obj <- CreateSeuratObject(counts = mtx,
                              assay = assay,
                              meta.data = adata$obs,
                              min.features = min_obs,
                              data = X,
                              min.cells = min_var,
                              project = project)
  s_obj[[assay]] <- AddMetaData(s_obj[[assay]], adata$var)
  
  # Save additional layers
  for(layer in names(adata$layers)){
    if(layer != counts){
      s_obj[[assay]][layer] <- dgR_flip(adata$layer[layer])
    }
  }
  
  # go through and add features. 
  if(!is.null(adata$var$highly_variable)){
    print("Adding Variable Features")
    VariableFeatures(s_obj) <- colnames(adata)[adata$var$highly_variable]
  }
  # Add PCA
  if(!is.null(adata$obsm$X_pca)){
    print("Adding Principal Components")
    pca_cells <- adata$obsm$X_pca
    rownames(pca_cells) <- rownames(adata)
    pca_genes <- adata$varm$PCs
    rownames(pca_genes) <- colnames(adata)
    pca <- CreateDimReducObject(
      embeddings = pca_cells,
      loadings = adata$varm$PCs,
      key = "PC_",
      assay = assay
    )
    s_obj[["pca"]] <- pca
  }
  # Add any embeddings
  embeddings <- c("TSNE", "UMAP")
  for(emb in embeddings){
    X_emb <- paste0("X_", tolower(emb))
    if(X_emb %in% names(adata$obsm)){
      print(paste("adding", emb, "embedding"))
      mtx <- adata$obsm[[X_emb]]
      rownames(mtx) <- rownames(adata)
      cell_embedding <- CreateDimReducObject(
        embeddings = mtx,
        key = paste0(emb, "_"),
        assay = assay
      )
      s_obj[[emb]] <- cell_embedding
    }
  }
  
  Idents(s_obj) <- idents
  
  if(export){
    saveRDS(s_obj, file=paste0(project, ".rds"))
  }
  
  return(s_obj)
}

Seurat2Scanpy <- function(seurat_object,
                          assay = "RNA",
                          X = "data",
                          backing_path = "adata.h5ad"){
  #' Convert SeuratObjects to Scanpy AnnData
  #' Retain as much information (ex: HVGs, embeddings) as possible.
  #' 
  #' @param seurat_object path to seurat_object RDS file, or SeuratObject
  #' @param assay assay to convert to AnnData
  #' @param X layer to use as main counts (adata$X)
  #' @param backing_path name of backing file
  
  # Read Seurat Object in
  if(class(s_obj) == "SeuratObject"){
    s_obj <- seurat_object
  } else if(is.character(seurat_object)){
    if(tools::file_ext(seurat_object) == "RDS"){
      s_obj <- readRDS(seurat_object)
    }
  } else {
    print(paste("Unable to read Seurat Object from", seurat_object))
    return(NA)
  }
  
  projections <- names(s_obj)[!(names(s_obj) %in% Assays(s_obj))]
  layers <- Layers(s_obj)[Layers(s_obj) != X]
  
  adata <- AnnData(
    X = NULL,
    obs = s_obj[[]],
    var = s_obj[[assay]],
    uns = NULL,
    obsm = lapply(projections, function(x){Embeddings(s_obj, x)}),
    varm = lapply(projections, function(x){Loadings(s_obj[x])}),
    layers = lapply(layers, function(x){s_obj[[assay]][x]}),
    raw = NULL,
    dtype = "float32",
    shape = NULL,
    filename = backing_path,
    filemode = NULL,
    obsp = NULL,
    varp = NULL
  )

  
  return(adata)
}

Scanpy2Loupe <- function(adata,
                         counts="X",
                         subset_column = NA,
                         clusters = c("leiden", "annot"),
                         projections = c("X_umap", "X_pca"),
                         out_dir = ".",
                         name = "scanpy"){
  #' @param adata AnnData object, or path to itdescription
  #' @param counts Layer containing counts, to be passed to Loupe Browser
  #' @param subset_column Break object across this column in .obs
  #' @param clusters Groups in .obs to pass as colors to Loupe Browser
  #' @param projections Projections in .obsm to pass to Loupe Browser
  #' @param out_dir Where to save cloupe object
  #' @param name Name cloupe will be saved under
  
  # Read Data
  if (is.character(adata)) {
    adata <- read_h5ad(adata)
  }
  
  if(is.na(subset_column)){
    mtx <- GetScanpyMtx(adata, counts)
    create_loupe(
      mtx,
      clusters = adata$obs[clusters],
      projections = lapply(adata$obsm[projections], function(x) x[,1:2]),
      output_dir = out_dir,
      output_name = name,
      feature_ids = rownames(adata$var),
      executable_path = NULL,
      force = FALSE,
      seurat_obj_version = "scanpy"
    )
  } else {
    subset_level = levels(adata$obs[, subset_column])
    
    for(lvl in subset_level){
      ad1 <- adata[adata$obs[, subset_column] == lvl]
      Scanpy2Loupe(ad1, 
                   counts=counts,
                   clusters=clusters, 
                   projections=projections,
                   out_dir = out_dir,
                   name=paste(name, lvl, sep="_"))
    }
  }
}



# ---- Command Line ----