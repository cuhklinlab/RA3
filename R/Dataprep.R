#' Normalize count matrix
#'
#' This funciton is optional for normalization of runRA3's input count matrix. It implements TF-IDF for matrix normailization.
#'
#'
#' @param sc_data Input single cell ATAC sequence count matrix, rows refer to features/regions, columns refer to cells.
#' @param ref_data Input reference data, rows refer to observations, columns refer to features/regions.
#' @return A list containing following components:
#' \item{sc_data}{Normalized scATAC-seq count matrix with selected features.}
#' \item{ref_data}{Normalized reference count matrix with selected features.}
#' @import pracma
#' @export
Dataprep <-  function(sc_data, ref_data){
  Y <- sc_data # p by n
  bulk_mat <- ref_data # n_bulk by p
  data_full <- list()

  # Data Preprocessing
  # TF-IDF
  nfreqs = Y / pracma::repmat(apply(Y, 2, sum), dim(Y)[1], 1)
  Y = nfreqs * t(pracma::repmat(log(1 + dim(Y)[2]) / apply(Y,1,sum), dim(Y)[2], 1))
  data_full$sc_data <- Y
  data_full$ref_data <- bulk_mat
  return(data_full)

}
