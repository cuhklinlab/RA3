# The R-based implementation of RA3
The R-based implementation of "A reference-guided approach for epigenetic characterization of single cells".<br/>
The source code for the reproduction of results in the manuscript can be found [here](https://github.com/cuhklinlab/RA3_source).

## Installation guide
Install the released version of **RA3** package using devtools from Github:
```r
install.packages("devtools")
devtools::install_github("cuhklinlab/RA3")
```

Package **Seurat** is needed for current version of **RA3**. Install **Seurat** using the following chunk, or check [here](https://satijalab.org/seurat/install.html).
```r
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
```

### Troubleshooting dependencies
At this point, there might be some missing dependicies from CRAN or Bioconductor. 

To install missing dependencies from Bioconductor, use the following chunk:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('chromVAR', 'motifmatchr', 'SummarizedExperiment', 'BSgenome.Hsapiens.UCSC.hg19', 'slingshot')) 
```

It will take a few minutes to install the **RA3** package, mainly for preparing the embedded demo-data. 

## Functions
This package includes following main functions:
- `runRA3` runs RA3 for the analysis of single-cell chromatin accessibility sequencing data. 
- `Dataprep` normalizes the input count matrix by TF-IDF.
- `RA3_EM` performs an EM algorithm to estimate parameters of the RA3 model.
- `RA3_clustering` deploys louvain clustering on the output of RA3.
- `RA3_TrajInf` does trajectory inference based on the output of RA3.
- `RA3_motif` runs motif enrichment based on the output of RA3.
- `RA3_pcatsne` deploys PCA and t-SNE for dimension reduction on normalized scCAS data.
- `RA3_RefProj` deploys a reference projection approach on the TF-IDF normalized scCAS data, then uses t-SNE for further dimension reduction.

## Documentation
Please check the [vigenette](https://github.com/cuhklinlab/RA3/wiki) for a tutorial. Two examples are contained for a quick start of RA3.

## System requirements
### Software requirements
The package has been tested on the following operating systems: 
- Linux: CentOS Linux release 7.7.1908
- macOS: Catalina (10.15.4)
- Windows: Windows 10

### Hardware requirements
The package has been tested on both Normal Personal Computer and High-Performance Computing Cluster.

## License
This package is built under license **GNU GENERAL PUBLIC LICENSE (GPL)**.
