#' Trajectory inference based on the latent features output by RA3
#'
#' The output of RA3 can also be implemented in trajectory inference. This function uses the extracted features from RA3 for trajectory inference by slingshot. 
#'
#' @param DRres_mat the extracted features H output by runRA3, colomns referring to cells
#' @param cluster_assign the estimated cluster labels
#' @param rand_seed random seed set for tsne plot, default value is set as 23
#' @return A list  containing following items:
#' \item{tsne}{matrix containing tsne representation, calculated using random seed 'rand_seed'}
#' \item{sds.new}{A SlingshotDataSet containing curves. This contains the curves representing the estimated cell lineages.}
#' @examples
#' result <- RA3_TrajInf(DRres_mat,cluster_assign,rand_seed = 23)
#' @importFrom Rtsne Rtsne
#' @importFrom slingshot getLineages getCurves embedCurves
#' @export
RA3_TrajInf <- function(DRres_mat,cluster_assign,rand_seed = 23){

  DRres_mat <- t(DRres_mat)
  result <- list()
  set.seed(rand_seed)
  result$tsne  <- Rtsne::Rtsne( DRres_mat )$Y 
  lin1 <- slingshot::getLineages(DRres_mat, cluster_assign)
  crv1 <- slingshot::getCurves(lin1)
  result$sds.new <- slingshot::embedCurves(crv1, result$tsne, shrink.method='density')

  return(result)

}
