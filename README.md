# RA3
RA3 method for integrative analysis of single cell ATAC sequence data.

## Requirements

## Installation

## Tutorial

This package includes following functions:
- `runRA3` Run RA3 for integrative analysis of scATAC-seq data. 
- `Dataprep` Normalize count matrix. This funciton is optional for normalization of runRA3's input count matrix. It consists of two steps:
  1. Use TF-IDF for matrix normailization,
  2. Select peaks for rough dimension reduction.
- `RA3_EM` EM algorithm for RA3. This function is implementing an EM algorithm to estimate parameters of RA3 model.
