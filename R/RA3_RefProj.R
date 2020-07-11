#' Dimension reduction using reference projection approach
#'
#' This function deploys a reference projection approach on the TF-IDF normalized scCAS data, then uses t-SNE for further dimension reduction.
#'
#' @param sc_data_nml TF-IDF normalized scCAS data, the output of function Dataprep
#' @return \item{bulk_tsne}{tsne representation for projected scCAS data}
#' @examples
#' result <- RA3_RefProj(ref_data_nml,sc_data_nml,rand_seed = 2020)
#' @importFrom Rtsne Rtsne
#' @export
RA3_RefProj <- function(ref_data_nml,sc_data_nml,rand_seed = 23){
  res_bulk <- prcomp(ref_data_nml, center = T, retx = T)
  coeff_bulk <- res_bulk$rotation
  set.seed(rand_seed)
  bulk_tsne <- Rtsne::Rtsne(t(sc_data_nml) %*% coeff_bulk)$Y
  bulk_tsne <- as.data.frame(bulk_tsne)
  colnames(bulk_tsne) <- c('tsne1','tsne2')
  return(bulk_tsne)
}