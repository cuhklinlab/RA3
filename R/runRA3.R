#' Run RA3 for integrative analysis of scATAC-seq data
#'
#' RA3 refers to Reference-guided Approach for the Analysis of scATAC-seq data. It can simultaneously incorporate shared biological variation from reference data and identify distinct subpopulations, and thus achieves superior performance to existing methods in comprehensive experiments.
#'
#' @param sc_data scATAC-seq count matrix, the rows should refer to features/regions and columns refer to cells.
#' @param ref_data reference data matrix, the columns should refer to features/regioins and rows refer to observations.
#' @param K2 the number of components in RA3's second part, the default value is K2 = 5.
#' @param K3 the number of components in RA3's third part, the default value is K3 = 5.
#' @param normalize a logical value indicating whether the output H should be normalized, the default value is set TRUE.
#' @param ttest a logical value indicating whether the output H2 should take a one sample t-test to select most informative components, the default value is set TRUE.
#' @return A list containing the following components:
#' \item{H}{the extracted latent features H.}
#' \item{W}{the estimated matrix of parameter matrix W.}
#' \item{Beta}{the estimated covariance parameter vector \eqn{\beta}.}
#' \item{Gamma}{the estimated indicator matrix \eqn{\Gamma}.}
#' \item{A}{the estimated precision matrix A.}
#' \item{sigma_s}{the estimated \eqn{\sigma^2}.}
#' \item{lgp}{the largest log posterior value when EM algorithm converges.}
#' @examples
#' result <- runRA3(sc_example, reference_example)
#' result <- runRA3(sc_example, reference_example, 10, 5)
#' @import pracma
#' @import irlba
#' @export

runRA3 <- function(sc_data, ref_data, K2 = 5, K3 = 5, normalize=TRUE, ttest=TRUE){
  Y <- sc_data # p by n
  bulk_mat <- ref_data # n_bulk by p

  # Data Preprocessing
  # peak selection with ratio 0.03
  filter_peak = apply((t(Y)>=1),2,sum) >= floor(dim(Y)[2]*0.03)
  Y = Y[filter_peak, ]
  bulk_mat = bulk_mat[ ,filter_peak]

  # TF-IDF
  nfreqs = Y / pracma::repmat(apply(Y, 2, sum), dim(Y)[1], 1)
  Y_mat = nfreqs * t(pracma::repmat(log(1 + dim(Y)[2]) / apply(Y,1,sum), dim(Y)[2], 1))

  # Calculate Initialization Value for the Algorithm
  # pca for reference
  pca_bulk <- prcomp(bulk_mat, center = T, retx = T)
  coeff <- pca_bulk$rotation
  score <- pca_bulk$x
  # Visualization of Reference Projection
  # bulk_tsne <- Rtsne(t(t(coeff) %*% Y_mat) )

  # Standardize
  p <- nrow(Y)
  n <- ncol(Y)
  K1 <- nrow(bulk_mat) - 1
  # K2 <- 5
  # K3 <- 5

  latent_h <- t(coeff) %*% Y_mat
  V_beta <- sqrt(apply(t(latent_h),2,var))
  init_V <- pracma::repmat(V_beta, p, 1) * coeff

  # Good Start of K2 Component
  residual <- Y_mat - coeff %*% (t(coeff) %*% Y_mat)
  # res_pca <- prcomp(t(residual), center = T, retx = T) # !!! number of components
  res_pca <- irlba::irlba(t(residual), nv = K2)
  coeff_res <- res_pca$v
  score_res <- res_pca$u %*% diag(res_pca$d)
  score_res_rotate <- varimax(score_res[ ,1:K2])
  RM <- score_res_rotate$rotmat
  coeff_res_rotate <- coeff_res[ ,1:K2] %*% RM # p by k
  score_reconst_rotate <- t(residual) %*% coeff_res_rotate # n by k
  # H2_ini <- score_reconst_rotate %*% diag(sqrt(apply(score_reconst_rotate, 2, var))^(-1))
  W2_ini <- coeff_res_rotate %*% diag(sqrt(apply(score_reconst_rotate, 2, var)))

  # Good Start of Sigma
  center_Y_stand <- Y_mat - t(rep(1, ncol(Y_mat)) %*% t(rowMeans(Y_mat)))
  residual_stand = center_Y_stand - coeff %*% (t(coeff) %*% center_Y_stand)
  # res_stand_pca <- prcomp(t(residual_stand), center = T, retx = T) # !!! number of components
  res_stand_pca <- irlba::irlba(t(residual_stand), nv = 20)
  coeff_res_stand <- res_stand_pca$v
  score_res_stand <- res_stand_pca$u %*% diag(res_stand_pca$d)
  score_rotate2 <- varimax(score_res_stand[ ,1:20])
  RM_stand <-  score_rotate2$rotmat
  coeff_stand_rotate <- coeff_res_stand[, 1:20] %*% RM_stand
  score_reconst_stand_rotate = t(residual_stand) %*% coeff_stand_rotate

  res_final <- residual_stand - coeff_stand_rotate %*% t(score_reconst_stand_rotate)
  epsilon_stard = res_final - t(rep(1, ncol(res_final)) %*% t(rowMeans(res_final)))
  # sigma_setting = psych::tr(t(epsilon_stard) %*% epsilon_stard)/(n*p)
  sigma_setting = sum(epsilon_stard^2)/(n*p)

  # Parameters Setting
  K <- K1 + K2 + K3
  para_num <- p * K


  # Initialization
  W1_ini <- init_V[ ,1:K1]
  # H1_ini <- t(latent_h[1:K1,]) %*% diag(V_beta[1:K1]^(-1))
  # H_PCA <- t(cbind(H1_ini, H2_ini[,1:K2], matrix(rnorm(n*K3,0,1), n, K3)))
  W_PCA <- cbind(W1_ini, W2_ini[ ,1:K2], matrix(rnorm(p*K3,0,1), p, K3))
  A_PCA <- diag(c(rep(1,K2), rep(1,K3)))
  Gamma_PCA <- matrix(rep(1,K*n), K, n)
  Gamma_PCA[(K1+1):(K2+K1), ] <-  matrix(rep(0, K2*n), K2, n)
 # V_PCA = cbind(init_V[ ,1:K1], matrix(rep(0,p*(K2+K3)), p, K2+K3))

  # Run
  # result <- RA3_EM(Y_mat,theta,sigma,sigma1,sigma2,K1,K2,K3,Gamma_PCA,A_PCA,W_PCA,sigma_setting,H_PCA,X,Beta,res_path)
   result <- RA3_EM(Y_mat,normalize,ttest,K1,K2,K3,Gamma_PCA,A_PCA,W_PCA,sigma_setting)
   return(result)
  }
