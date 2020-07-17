#' Cell clustering based on the latent features output by RA3
#'
#' Implement louvain clustering based on the output of RA3, using package Seurat.
#'
#' @param H the extracted features H output by runRA3, colomns referring to cells
#' @param n_cluster the number of clusters
#' @return \item{Seurat_louvain}{vector containing cell cluster assignment}
#' @examples
#' result <- RA3_clustering(H, n_cluster)
#' @export
RA3_clustering <- function(H, n_cluster){
  rownames(H) <- as.character(1:nrow(H))
  colnames(H) <- as.character(1:ncol(H))
  H_obj <- Seurat::CreateSeuratObject(counts = H)
  H_obj <- Seurat::ScaleData(H_obj, features = rownames(H_obj))
  H_obj <- Seurat::RunPCA(H_obj, features = rownames(H_obj))
  H_SNN <- Seurat::FindNeighbors(H_obj, k.param = 20, features = rownames(H_obj), dims = 1:(nrow(H)-1))
  # clusterlouvain2 <- getNClusters(H_SNN, n_clusters = n_cluster) 
  
  min_res <- 0
  max_res <- 3
  max_steps <- 20
  this_step <- 0
  while (this_step < max_steps) {
    cat('step', this_step,'\n')
    this_res <- min_res + ((max_res-min_res)/2)
    data_louvain <- Seurat::FindClusters(H_SNN, resolution = this_res, verbose = FALSE)
    this_clusters <- length(unique(Seurat::Idents(data_louvain)))
    
    if (this_clusters > n_cluster) {
      max_res <- this_res
    } else if (this_clusters < n_cluster) {
      min_res <- this_res
    } else {
      break
    }
    
    this_step <- this_step + 1    
  }
  
  
  Seurat_louvain <- Seurat::Idents(data_louvain)
  return(Seurat_louvain)
}



