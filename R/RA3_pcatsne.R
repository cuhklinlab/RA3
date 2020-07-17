#' Dimension reduction using PCA and t-SNE
#'
#' Implement PCA on the TF-IDF normalized scCAS data, then deploy t-SNE on the component scores.
#'
#' @param sc_data_nml TF-IDF normalized scCAS data, the output of function Dataprep
#' @return \item{pca_tsne}{tsne representation for component scores of the TF-IDF normalized scCAS data}
#' @examples
#' result <- RA3_pcatsne(sc_data_nml)
#' @importFrom irlba irlba
#' @importFrom Rtsne Rtsne
#' @export
RA3_pcatsne <- function(sc_data_nml, rand_seed = 23){
  res_pca <- irlba::irlba(t(sc_data_nml), nv=10)
  score_pca <- res_pca$u %*% diag(res_pca$d)
  set.seed(rand_seed)
  pca_tsne <- Rtsne::Rtsne(score_pca)$Y
  pca_tsne <- as.data.frame(pca_tsne)
  colnames(pca_tsne) <- c('tsne1','tsne2')
  return(pca_tsne)
}