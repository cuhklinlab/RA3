# RA3
RA3 method for integrative analysis of single cell ATAC sequence data.

## Installation
You can install the released version of package Ra3 from Github:
```javascript
devtools::install_github("cuhklinlab/RA3")
```
## Main funcitons

This package includes following functions:
- `runRA3` Run RA3 for integrative analysis of scATAC-seq data. 
- `Dataprep` Normalize count matrix. This funciton is optional for normalization of runRA3's input count matrix. It consists of two steps:
  1. Use TF-IDF for matrix normailization,
  2. Select peaks for rough dimension reduction.
- `RA3_EM` EM algorithm for RA3. This function is implementing an EM algorithm to estimate parameters of RA3 model.

## Tutorial
Please check the [vigenette](https://github.com/cuhklinlab/RA3/wiki) for a tutorial. Two examples are contained for a quick start of RA3.
