#' EM algorithm for RA3
#'
#' This function is implementing an EM algorithm to estimate parameters of RA3 model.
#'
#' @param Y Input single cell count matrix.
#' @param normalize A logical value indicating whether the output W should be normalized through columns. Corresponding scale would be multiplied onto H if \emph{normalize} is TRUE.
#' @param ttest A logical value indicating whether a one sample ttest should be done for output H2 components.
#' @param K1,K2,K3 Number of components consisted in RA3 model.
#' @param Gamma Initial matrix of parameter matrix \eqn{\Gamma}.
#' @param A Initial matrix of precision matrix A.
#' @param W Initial matrix of parameter matrix W, a warm start is recommended in \code{\link{runRA3}}.
#' @param sigma_s Intial value of parameter \eqn{\sigma^2}, a warm start is recommended in \code{\link{runRA3}}.
#' @return A list containing the following components:
#' \item{H}{the estimated latent variable H.}
#' \item{W}{the estimated matrix of parameter matrix W.}
#' \item{Beta}{the estimated covariance parameter vector \eqn{\beta}.}
#' \item{Gamma}{the estimated indicator matrix \eqn{\Gamma}.}
#' \item{A}{the estimated precision matrix A.}
#' \item{sigma_s}{the estimated \eqn{\sigma^2}.}
#' \item{lgp}{the largest log posterior value when EM algorithm converges.}
#'
#' @import pracma
RA3_EM <- function(Y,normalize,ttest,K1,K2,K3,Gamma,A,W,sigma_s){
  pmt <- proc.time()[3]
  res <- list()

  # setting
  p <- nrow(Y)
  n <- ncol(Y)
  K <- K1 + K2 + K3
  theta <- 0.1
  tau <- 1
  tau0 <- 0.9
  tau1 <- 5
  q <- 1
  X <- matrix(1,1,n)
  Beta <- matrix(0,p,q)
  if (is.matrix(sigma_s)){
    sigma_s <- sigma_s[1]
  }
  MAXIte <- 5000
  err <- 1e-5

  Qy <- rep(-Inf, MAXIte)
  Qw <- rep(-Inf, MAXIte)
  Qh <- rep(-Inf, MAXIte)
  Qgamma <- rep(-Inf, MAXIte)
  log_p <- rep(-Inf, MAXIte)
  diff_log_Q <- rep(-Inf, MAXIte)
  Q_H <- rep(-Inf, MAXIte)
  max_A <- c(1e12*rep(1,K1), 1e12*rep(1,K2), 1e12*rep(1,K3))
  cont <- rep(0, n)
  phi_gamma_0 <- tau0^-2
  phi_gamma_1 <- tau1^-2

  # Run
  for(ite in 2:MAXIte){
    # E step
    E_h <- matrix(rep(0, K*n),K,n)
    E_hh <- array(rep(0, K*K*n), c(K,K,n))
    W_s <- 1/sigma_s * t(W) %*% W
    D <- (matrix(rep(1, K*n),K,n) - Gamma) * tau0^2 + Gamma * tau1^2
    D[1:K1, ] <- tau^2 * matrix(rep(1, K1*n),K1,n)
    D[(K1+K2+1):K, ] <- tau^2 * matrix(rep(1, K3*n),K3,n)
    for (j in 1:n){
        D_j <- D[ ,j]
        SIGMA_h_inv <- W_s + diag(D_j^(-1)) # K*K
        SIGMA_h <- solve(SIGMA_h_inv) # K*K
        MU_h <- 1/sigma_s *t(Y[,j]-Beta %*% X[,j]) %*% W %*% SIGMA_h
        E_h[,j] <- t(MU_h)
        E_hh[ , ,j]  <- E_h[,j] %*% t(E_h[,j]) + SIGMA_h
    }

    # M step
    W[ ,(K1+1) : K] <- (Y - Beta %*% X - W[ ,(1:K1)] %*% E_h[(1:K1), ]) %*% t(E_h[(K1+1) : K, ]) %*% solve(apply(E_hh[(K1+1):K, (K1+1):K, ],c(1,2),sum) + sigma_s * A, tol = 1e-100)
    Beta <- (Y - W %*% E_h) %*% t(X) %*% solve(X %*% t(X))
    A_new <- diag(A)
    for (k in 1:(K2+K3)){
      tmp <- W[ ,k+K1]
      A_new[k] <- p / (t(tmp) %*% tmp)
      if (A_new[k] > max_A[k]){
              A_new[k] <- max_A[k]
      }
    }
    A <- diag(A_new)

    Gamma_new <- matrix(rep(1,K*n),K,n)
    for (k in (K1+1):(K1+K2)){
    for (j in 1:n){
      f_gamma_0 <- - 0.5 * E_hh[k,k,j] * phi_gamma_0 + 0.5 * log(phi_gamma_0) + log(1-theta)
      f_gamma_1 <- - 0.5 * E_hh[k,k,j] * phi_gamma_1 + 0.5 * log(phi_gamma_1) + log(theta)
      Gamma_new[k,j] <- 1 * (f_gamma_1 > f_gamma_0)
      }
     }
    Gamma <- Gamma_new

    sigma_s <- 0
    W_tmp <- t(W) %*% W
    TMP <- Y - Beta %*% X
    Mu_tmp <- t(W) %*% TMP
    for (j in 1:n){
      # sigma_s <- sigma_s + t(TMP[,j]) %*% TMP[,j] - 2*t(Mu_tmp[,j]) %*% E_h[,j] + psych::tr(E_hh[ , ,j] %*% W_tmp)
      sigma_s <- sigma_s + t(TMP[,j]) %*% TMP[,j] - 2*t(Mu_tmp[,j]) %*% E_h[,j] + sum(E_hh[ , ,j] * t(W_tmp))
    }
    sigma_s <- sigma_s[1] / (p*n)

    # Calculate log posterior
    # Qw[ite] <- - (K2+K3) * p/2 * log(2*pi) + p/2 * psych::tr(log(A)) - 1/2 * psych::tr(A %*% t(W[ ,(K1+1):K]) %*% W[ ,(K1+1):K])
    Qw[ite] <- - (K2+K3) * p/2 * log(2*pi) + p/2 * sum(diag(log(A))) - 1/2 * sum(A %*% t(W[ ,(K1+1):K]) * t(W[ ,(K1+1):K]))

    Qgamma[ite] <- sum(Gamma[(K1+1):(K1+K2), ] * log(theta) + (matrix(rep(1,K2*n),K2,n) - Gamma[(K1+1):(K1+K2),]) * log(1-theta))
    log_tmp <- - n * p/2 * log(2*pi*sigma_s)
    W_s <- 1/sigma_s * W_tmp
    D <- (matrix(rep(1,K*n),K,n) - Gamma) * tau0^2 + Gamma * tau1^2
    D[1:K1, ] <- tau * matrix(rep(1,K1*n),K1,n)
    D[(K1+K2+1):K, ] <- tau * matrix(rep(1,K3*n),K3,n)

    for (j in 1:n){
      D_j <- D[ ,j]
      Sigma_j_inv <- W_s + diag(D_j^(-1))
      Sigma_j <- solve(Sigma_j_inv)
      Mu_j <- 1/sigma_s *  Sigma_j %*% Mu_tmp[ ,j]
      # cont[j] <- - 1/2 * log(prod(D_j)) - 1/2 * (1/sigma_s * t(TMP[ ,j]) %*% TMP[ ,j] ) + 1/2 * log(det(Sigma_j)) + 1/2 * t(Mu_j) %*% solve(Sigma_j) %*% Mu_j
      cont[j] <- - 1/2 * log(prod(D_j)) - 1/2 * 1/sigma_s * t(TMP[,j]) %*% TMP[,j] + 1/2 * log(det(Sigma_j)) + 1/2 * t(Mu_j) %*% Sigma_j_inv %*% Mu_j
      log_tmp <- log_tmp + cont[j]
    }
    Qh[ite] <- log_tmp
    log_p[ite] <- log_tmp + Qw[ite] + Qgamma[ite]

    if(log_p[ite] > log_p[ite-1]){
      res$H <- E_h
      res$W <- W
      res$Beta <- Beta
      res$Gamma <- Gamma
      res$A <- A
      res$sigma_s <- sigma_s
      res$lgp <- log_p[ite]
      } else
      stop('[Warning] Reduced log_p.',call = FALSE)

    pracma::fprintf('[%d] Time: %0.1fs\tH: %s\tW: %s\tGamma: %s\tlog_p: %s\n',ite,proc.time()[3]-pmt,as.character(Qh[ite]),as.character(Qw[ite]),as.character(Qgamma[ite]),as.character(log_p[ite]), file = "" )

    if(abs(log_p[ite] - log_p[ite-1]) < err * abs(log_p[ite-1])){
      break
    }
  }

  # Normalize
  if (normalize == TRUE) {
    scale <- sqrt(diag(t(res$W) %*% res$W))
    res$W <- res$W / pracma::repmat(scale,p,1)
    res$H <- res$H * matrix(scale,K,n)
  }

  # t_test for outputing H
  if (ttest == TRUE) {
    H2_ind <- rep(0,K2)
    for (k in 1:K2) {
      H_ttest <- stats::t.test(res$H[(K1+k), ])
      if (H_ttest$p.value <= 0.05){
        H2_ind[k] <- 1
      } else
        H2_ind[k] <- 0
    }
    H_ind <- c(rep(1,K1), H2_ind, rep(1,K3))
    res$H <- res$H[which(H_ind == 1), ]
  }

  pracma::fprintf('\nSetting\nN: %d\tP: %d\tK: %d\tTime: %0.1fs\tmax_log_p: %s\n',n,p,K,proc.time()[3]-pmt,as.character(max(log_p)), file = "")
  return(res)
}

# Other Functions
# norm_vec <- function(x) sqrt(sum(x^2))
