library(reticulate)
library(anndata)
library(Seurat)
library(Matrix)

py_require(c("anndata"))

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

Scanpy2Seurat <- function(h5ad,
                          counts = "X",
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
  
  if(is.environment(h5ad)){
    adata <- h5ad
  } else if (is.character(h5ad)) {
    adata <- read_h5ad(h5ad)
  } else {
    print("Unable to read h5ad file. Please try again.")
    quit()
  }
  
  # Get counts and transform to match SeuratObject conventions
  if( is.null(counts) || (counts == "X") ){
    mtx <- dgR_flip(adata$X)
    X <- NULL
  } else {
    mtx <- dgR_flip(adata$layers[counts])
    X <- dgR_flip(adata$X)
  }
  
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
