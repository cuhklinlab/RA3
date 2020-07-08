# The R-based implementation of RA3
The R-based implementation of "A reference-guided approach for epigenetic characterization of single cells".<br/>
The source code for the reproduction of results in the manuscript can be found [here](https://github.com/cuhklinlab/RA3_source).

## Installation guide
Install the released version of RA3 package from Github:
```javascript
devtools::install_github("cuhklinlab/RA3")
```

It will take a few minutes to install the **RA3** package, mainly for preparing the embedded demo-data. 

## Funcitons
This package includes following main functions:
- `runRA3` runs RA3 for the analysis of single-cell chromatin accessibility sequencing data. 
- `Dataprep` normalizes the input count matrix by TF-IDF.
- `RA3_EM` performs an EM algorithm to estimate parameters of the RA3 model.

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
