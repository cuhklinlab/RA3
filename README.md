# RA3
RA3 method for integrative analysis of single cell ATAC sequence data.

## System requirements
### Hardware requirements
'RA3' package only needs a standard computer with enough RAM to support the operations defined by a user.

### Software requirements
The package has been tested on the following systems: 
- Linux:
- macOS: Catalina (10.15.4), 
- Windows: Windows 10

## Installation guide
You can install the released version of package RA3 from Github:
```javascript
devtools::install_github("cuhklinlab/RA3")
```

It will take around 1 minute to install package RA3 on a standard MacBook Pro computer (8 GB RAM, 4 cores@2.4 GHz) and internet of speed 50 Mbps. 

## Funcitons
This package includes following functions:
- `runRA3` Run RA3 for integrative analysis of scATAC-seq data. 
- `Dataprep` Normalize count matrix. This funciton is optional for normalization of runRA3's input count matrix. It consists of two steps:
  1. Use TF-IDF for matrix normailization,
  2. Select peaks for rough dimension reduction.
- `RA3_EM` EM algorithm for RA3. This function is implementing an EM algorithm to estimate parameters of RA3 model.

## Documentation
Please check the [vigenette](https://github.com/cuhklinlab/RA3/wiki) for a tutorial. Two examples are contained for a quick start of RA3.

## License
This package is built under license **GNU GENERAL PUBLIC LICENSE (GPL)**.
