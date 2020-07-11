#' getbetaMLE
#'
#' @param data
#' @param cluster
#' @param offset
#' @param beta_ini
#' @param maxiter
#' @param thres
#' @param quiet
#' @return beta
getbetaMLE <- function(data, cluster, offset, beta_ini, maxiter, thres, quiet){
  ## data is n by p
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## offset corresponds to h_i, here it is a vector of length n
  p <- ncol(data)
  x <- getdesign(cluster)
  beta <- beta_ini
  betaDiffs <- c() 
  percConverged <- c()
  converged_flag <- rep(0, p)
  
  converged <- 0
  iter <- 1
  while (!converged & iter<maxiter){
    mu <- offset*exp(x%*%beta) 
    u <- crossprod(x, data-mu)
    updateBeta <- function(r){
      if (converged_flag[r]==0 & !is.na(beta[1,r]) ){
        Jr <- crossprod(x, mu[,r]*x)
        betaDiff <- try(solve(Jr, u[,r]), silent=T)
        if (class(betaDiff)=="try-error"){
          return( rep(NA, length(beta[,r])) )    
        } else {
          return(beta[,r] + betaDiff)  
        }
      } else {
        return(beta[,r])
      }
    }
    
    betaPre <- beta
    beta <- sapply(1:p, updateBeta)
    betaDiffs <- rbind(betaDiffs, colSums(abs(beta - betaPre)) )
    converged_flag <- (betaDiffs[iter,] <= thres) + 0
    converged_flag[which(is.na(converged_flag))] <- 0
    singular_flag <- is.na(beta[1,]) + 0
    percConverged <- c(percConverged, sum(converged_flag[which(singular_flag==0)])/sum(singular_flag==0) )
    
    if (!quiet){
      cat("\r", round(percConverged[length(percConverged)]*100), "%", "converged")  
    }
    if (percConverged[length(percConverged)]==1){
      converged <- 1
    }
    iter <- iter + 1
  }
  return(beta)
}

#' getSigmaPrior
#'
#' @param beta
#' @param q
#' @return x
getSigmaPrior <- function(beta, q=0.05){
  ## get the empirical prior estimate
  ## q is the quantile to match
  
  ## center each column in beta
  betaC <- t(t(beta) - colMeans(beta))
  ## 
  sigmas <- rep(0, nrow(betaC))
  for (i in 1:nrow(betaC)){
    tmp <- betaC[i,]
    tmp <- tmp[which(!is.na(tmp))]
    sigmas[i] <- quantile(abs(tmp), 1-q)/qnorm(1-q/2)
  }
  return(sigmas)
}

#' getdesign
#'
#' @param cluster
#' 
getdesign <- function(cluster){
  nCluster <- length(unique(cluster))
  n <- length(cluster)
  x <- matrix(0, nrow=n, ncol=nCluster)
  for (i in 1:nCluster){
    x[which(cluster==i), i] <- 1
  }
  return(x)
}

#' getContrast
#'
#' @param nCluster
#' @return A list containing following items:
#' \item{constrasts}{contrasts needed to calculate the p-value}
#' \item{constrastsCluster}{clusters for contrasts}
getContrast <- function(nCluster){
  ## get all the contrasts that is needed to calculate the p-value
  constrasts <- matrix(0, nrow=(nCluster-1)*nCluster, ncol=nCluster)
  constrastsCluster <- rep(1:nCluster, each=nCluster-1)
  for (i in 1:nCluster){
    constrasttmp <- matrix(0, nrow=nCluster-1, ncol=nCluster)
    constrasttmp[, i] <- 1
    count <- 1
    for (j in 1:(nCluster-1)){
      constrasttmp[count, (c(1:nCluster)[-i])[j]] <- -1
      count <- count + 1
    }
    constrasts[(i*(nCluster-1)-(nCluster-1)+1):(i*(nCluster-1)), ] <- constrasttmp
  }
  return(list(constrasts=constrasts, constrastsCluster=constrastsCluster))
}

#' getbetaMAP
#'
#' @param data
#' @param cluster
#' @param offset
#' @param sigma
#' @param beta_ini
#' @param maxiter
#' @param thres
#' @param quiet
#' @return A list of beta and p-value
getbetaMAP <- function(data, cluster, offset, sigmas, beta_ini, maxiter, thres, quiet){
  ## data is n by p
  ## beta_ini includes the intercept
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## offset corresponds to h_i, here it is a vector of length n
  
  p <- ncol(data)
  x <- getdesign(cluster)
  ## add the intercept in x
  x <- cbind(1, x)
  ## get lambda
  lambda <- c(0, 1/sigmas^2)
  
  beta <- beta_ini
  betaDiffs <- c()
  percConverged <- c()
  converged_flag <- rep(0, p)
  
  converged <- 0
  iter <- 1
  while (!converged & iter<maxiter){
    mu <- offset*exp(x%*%beta) 
    z <- log(mu/offset) + (data-mu)/mu
    
    updateBetaRidge <- function(r){
      if (converged_flag[r]==0){
        JrRidge <- crossprod(x, mu[,r]*x) + diag(lambda)
        betar <- solve(JrRidge, crossprod(x, matrix(mu[,r]*z[,r], ncol=1)) )
        return( as.vector(betar) )
      } else {
        return(beta[,r])
      }
    }
    
    betaPre <- beta
    beta <- sapply(1:p, updateBetaRidge)
    betaDiffs <- rbind(betaDiffs, colSums(abs(beta - betaPre)) )
    converged_flag <- (betaDiffs[iter,] <= thres) + 0
    percConverged <- c(percConverged, sum(converged_flag)/p )
    
    if (!quiet){
      cat("\r", round(percConverged[length(percConverged)]*100), "%", "converged")  
    }
    
    if (percConverged[length(percConverged)]==1){
      converged <- 1
    }
    iter <- iter + 1
  }
  mu <- offset*exp(x%*%beta) 
  
  ## get all the contrast matrices and cluster label for each contrast
  nCluster <- length(unique(cluster))
  tmp <- getContrast(nCluster)
  constrasts <- tmp$constrasts
  
  ## pad the intercept with 0 in the constrast
  constrasts <- cbind(0, constrasts)
  constrastsCluster <- tmp$constrastsCluster
  
  calRidgePvalueA <- function(r){
    ## calculates the p-value for a small hypothesis
    if ( !is.na(beta[1,r]) ){
      Jr <- crossprod(x, mu[,r]*x)
      JrRidgeInv <- solve(Jr + diag(lambda))
      covRidge <- JrRidgeInv%*%Jr%*%JrRidgeInv
      ##
      betaC <- constrasts%*%matrix(beta[,r], ncol=1)  
      CcovC <- constrasts%*%covRidge%*%t(constrasts)
      SEbetaC <- sqrt(diag(CcovC))
      ##
      pvs <- pnorm(betaC/SEbetaC, lower.tail = FALSE)
      return(pvs)
    } else {
      return(rep(NA, length(constrastsCluster)))
    }
  }
  if (!quiet){
    cat("\nCalculating the p-values\n")
  }
  pvsA <- sapply(1:p, calRidgePvalueA) 
  ## get the p-value for the large null hypothesis by taking maximum
  if (nCluster>2){
    pvs <- c()
    for (i in 1:nCluster){
      pvs <- cbind(pvs, apply(pvsA[which(constrastsCluster==i),], 2, max))
    }  
  } else {
    pvs <- t(pvsA)
  }
  return(list(beta=beta, pvalue=pvs))
}

#' getClusterSpecificPvalue
#'
#' @param data
#' @param cluster
#' @param offset
#' @param landmark
#' @param maxiter
#' @param thresMLE
#' @param thresMAP
#' @param quiet
#' @return A list of beta, p-value and sigmas.
getClusterSpecificPvalue <- function(data, cluster, offset, landmark=NULL, maxiter=1000, thresMLE=10^-3, thresMAP=10^-5, quiet=FALSE){
  ## the main function for peak selection
  ## data is nPeaks(p) by Cells(n)
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## offset corresponds to h_i, here it is a vector of length n
  ## landmark is optional. If landmark is provided, we only do hypothesis testing on the union of landmark peaks
  if (!is.null(landmark)){
    ## take union of the landmark peaks
    unionlandmark <- which(rowSums(landmark)!=0)
    p <- nrow(data)
    data <- data[unionlandmark,]
  } else {
  }
  
  data <- t(data) # make data n by p
  ## remove the offset=0 samples
  data <- data[which(offset>0),]
  cluster <- cluster[which(offset>0)]
  offset <- offset[which(offset>0)]
  ## get beta_MLE
  if (!quiet){
    cat("\nEstimating beta MLE\n")  
  }
  betaMLE_ini <- matrix(0, nrow=length(unique(cluster)), ncol=ncol(data))
  betaMLE <- getbetaMLE(data=data, cluster=cluster, offset=offset, beta_ini=betaMLE_ini, maxiter=maxiter, thres=thresMLE, quiet=quiet)
  
  ## get the empirical prior sigma
  sigmas <- getSigmaPrior(betaMLE)
  
  ## get beta_MAP and the p-value
  if (!quiet){
    cat("\nEstimating beta MAP\n")
  }
  betaMAP_ini <- rbind(0, matrix(0, nrow=length(unique(cluster)), ncol=ncol(data)))
  result <- getbetaMAP(data=data, cluster=cluster, offset=offset, sigmas=sigmas, beta_ini=betaMAP_ini, maxiter=maxiter, thres=thresMAP, quiet=quiet)
  
  if (!is.null(landmark)){
    betaMAP <- matrix( nrow=length(unique(cluster))+1, ncol=p )
    pvalue <- matrix( nrow=p, ncol=length(unique(cluster)) )
    betaMAP[, unionlandmark] <- result$beta
    pvalue[unionlandmark,] <- result$pvalue
  } else {
    betaMAP <- result$beta
    pvalue <- result$pvalue
  }
  colnames(pvalue) <- paste("Cluster", 1:length(unique(cluster)) )
  return(list(betaMAP=betaMAP, sigmas=sigmas, pvalue=pvalue))
}





#' Motif enrichment analysis based on the latent features output by RA3
#'
#' The output of RA3 can also be implemented in motif enrichment analysis. This function uses the extracted features from RA3 for motif analysis though package chromVAR . 
#'
#' @param peaks a vector containing peak information of the scCAS data, with elements formatted as 'chr1_24516561_24517061'
#' @param label_true true cell labels
#' @param label_est the estimated cluster labels, starting from 1
#' @param sc_mat scCAS count matrix
#' @param cluster_peak_num number of peaks selected for each cluster, default value is 1000
#' @param motif_num number of motifs shown in the plot, default value is 50
#' @return A list containing following items:
#' \item{metaDataInd}{cell index grouped by clustering results and cell labels}
#' \item{top_devs}{derivation scores estimated using chromVAR}
#' @examples
#' result <- RA3_motif(peaks, label_true, label_est, sc_mat,cluster_peak_num=1000, motif_num=50)
#' 
#' @importFrom chromVAR addGCBias getJasparMotifs computeDeviations deviationScores computeVariability
#' @importFrom motifmatchr matchMotifs
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @export
RA3_motif <- function(peaks, label_true, label_est, sc_mat, cluster_peak_num=1000, motif_num=50){

  result<- list()
  
  metaData <- data.frame(true_label=(label_true), pred_label=label_est)
  metaDataInd = with(metaData, order(pred_label, true_label))

  resultPeakSelection <- getClusterSpecificPvalue(data=sc_mat, cluster=label_est, offset=colSums(sc_mat))
  pvaluematrix=resultPeakSelection$pvalue
  rownames(pvaluematrix)=peaks
  colData <- data.frame(true_label=(label_true))
  
  colData <- colData[metaDataInd,]
  
  peakselected = c()
  for (i in 1:ncol(pvaluematrix)) {
    peakselected = c(peakselected,head(order(pvaluematrix[,i], decreasing = FALSE), cluster_peak_num))
  }
  peakselected = unique(sort(peakselected))
  
  chr=c()
  p1=c()
  p2=c()
  for (i in 1:length(peaks[peakselected])) {
    strings=unlist(strsplit(peaks[peakselected][i], "_"))
    chr=c(chr,strings[1])
    p1=c(p1,strings[2])
    p2=c(p2,strings[3])
  }
  
  peakrange=data.frame(chr,p1,p2)
  peakrange$p1=as.numeric(as.character(peakrange$p1))
  peakrange$p2=as.numeric(as.character(peakrange$p2))
  peaks.gr = GenomicRanges::GRanges(peakrange[,1], IRanges::IRanges(peakrange[,2], peakrange[,3]))
  
  
  
  matrix_use=sc_mat
  matrix_select=matrix_use[peakselected,] 
  matrix_select = matrix_select[,metaDataInd]
  colnames(matrix_select)=1:ncol(matrix_select)
  
  
  frag_counts = SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=matrix_select), rowRanges=peaks.gr, colData=colData)
  frag_counts = chromVAR::addGCBias(frag_counts, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  motifs <- chromVAR::getJasparMotifs()
  motifs.matched = motifmatchr::matchMotifs(motifs, frag_counts, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  dev = chromVAR::computeDeviations(object = frag_counts, annotations = motifs.matched)
  dev.scores = chromVAR::deviationScores(dev)
  variability = chromVAR::computeVariability(dev)
  
  
  
  top_motifs = variability$name[head(order(variability$variability, decreasing = TRUE), motif_num)]
  names(top_motifs) = rownames(variability[head(order(variability$variability, decreasing = TRUE), motif_num), ])
  top_devs = dev.scores[which(rownames(dev.scores) %in% names(top_motifs)), ]
  rownames(top_devs) = top_motifs[match(rownames(top_devs), names(top_motifs))]
  
  result$metaDataInd <- metaDataInd
  result$top_devs <- top_devs
  return(result)
}


