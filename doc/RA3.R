## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(RA3)

## -----------------------------------------------------------------------------
library(RA3)

## -----------------------------------------------------------------------------
sc_data <- donorBM0828_sc_mat
ref_data <- donorBM0828_bulk_mat

## ----eval=TRUE, results='hide', message=FALSE---------------------------------
result <- runRA3(sc_data, ref_data)

## ---- eval=TRUE---------------------------------------------------------------
data_nml <- Dataprep(sc_data, ref_data)
sc_data_nml <- data_nml$sc_data
ref_data_nml <- data_nml$ref_data

## ---- eval=TRUE, fig.align='center', fig.height=4, fig.width=5----------------
library(ggplot2)
res_pca <- irlba::irlba(t(sc_data_nml), nv=10)
score_pca <- res_pca$u %*% diag(res_pca$d)
pca_tsne <- Rtsne::Rtsne(score_pca)$Y
pca_tsne <- as.data.frame(pca_tsne)
colnames(pca_tsne) <- c('tsne1','tsne2')
cell_label <- donorBM0828_label_mat
ggplot(pca_tsne) + geom_point(aes(tsne1, tsne2, colour = cell_label), shape = 10, position = ggplot2::position_jitter()) + labs(title = "pca tsne")


## ---- eval=TRUE, fig.align='center', fig.height=4, fig.width=5----------------
res_bulk <- prcomp(ref_data_nml, center = T, retx = T)
coeff_bulk <- res_bulk$rotation
bulk_tsne <- Rtsne::Rtsne(t(sc_data_nml) %*% coeff_bulk)$Y
bulk_tsne <- as.data.frame(bulk_tsne)
colnames(bulk_tsne) <- c('tsne1','tsne2')
ggplot(bulk_tsne) + geom_point(aes(tsne1, tsne2, colour = cell_label), shape = 10, position = ggplot2::position_jitter()) + labs(title="bulk projection tsne")

## ---- eval=TRUE, fig.align='center', fig.height=4, fig.width=5----------------
H_tsne <- Rtsne::Rtsne(t(result$H[1: (nrow(result$H)-5), ]))$Y
H_tsne <- as.data.frame(H_tsne)
colnames(H_tsne) <- c('tsne1', 'tsne2')
ggplot(H_tsne) + geom_point(aes(tsne1, tsne2, colour = cell_label), shape = 10, position = ggplot2::position_jitter()) + labs(title = "RA3 tsne")


## ---- eval=TRUE,message=FALSE,results='hide', fig.height=4, fig.width=5-------
sc_data2 <- forebrain_sc_mat
ref_data2 <- forebrain_bulk_mat

## ---- eval=TRUE, results='hide', message=FALSE, fig.height=4, fig.width=5-----
result2 <- runRA3(sc_data2, ref_data2)

## ---- eval=TRUE---------------------------------------------------------------
data2 <- Dataprep(sc_data2, ref_data2)
sc_data2_nml <- data2$sc_data
ref_data2_nml <- data2$ref_data

## ---- eval=TRUE, fig.align='center', fig.height=4, fig.width=5----------------
library(ggplot2)
res2_pca <- irlba::irlba(t(sc_data2_nml), nv=10)
score2_pca <- res2_pca$u %*% diag(res2_pca$d)
pca2_tsne <- Rtsne::Rtsne(score2_pca)$Y
pca2_tsne <- as.data.frame(pca2_tsne)
colnames(pca2_tsne) <- c('tsne1','tsne2')
cell_label2 <- forebrain_label_mat
ggplot(pca2_tsne) + geom_point(aes(tsne1, tsne2, colour = cell_label2), shape = 10, position = ggplot2::position_jitter()) + labs(title = "pca tsne")


## ---- eval=TRUE, fig.align='center', fig.height=4, fig.width=5----------------
res2_bulk <- prcomp(ref_data2_nml, center = T, retx = T)
coeff2_bulk <- res2_bulk$rotation
bulk2_tsne <- Rtsne::Rtsne(t(sc_data2_nml) %*% coeff2_bulk)$Y
bulk2_tsne <- as.data.frame(bulk2_tsne)
colnames(bulk2_tsne) <- c('tsne1','tsne2')
ggplot(bulk2_tsne) + geom_point(aes(tsne1, tsne2, colour = cell_label2), shape = 10, position = ggplot2::position_jitter()) + labs(title="bulk projection tsne")

## ---- eval=TRUE, fig.align='center', fig.height=4, fig.width=5----------------
H2_tsne <- Rtsne::Rtsne(t(result2$H[1: (nrow(result2$H)-5), ]))$Y
H2_tsne <- as.data.frame(H2_tsne)
colnames(H2_tsne) <- c('tsne1', 'tsne2')
ggplot(H2_tsne) + geom_point(aes(tsne1, tsne2, colour = cell_label2), shape = 10, position = ggplot2::position_jitter()) + labs(title = "RA3 tsne")

