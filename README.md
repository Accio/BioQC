[![Build Status](https://travis-ci.org/Accio/BioQC.svg?branch=master)](https://travis-ci.org/Accio/BioQC)
[![codecov](https://codecov.io/gh/Accio/BioQC/branch/master/graph/badge.svg)](https://codecov.io/gh/Accio/BioQC)


# BioQC: Detect tissue heterogeneity in gene expression data

BioQC is a R/[Bioconductor](https://bioconductor.org/packages/release/bioc/html/BioQC.html) package to detect tissue heterogeneity in gene expression data. To this end BioQC implements a computationally efficient Wilcoxon-Mann-Whitney test and provides more than 150 signatures of tissue-enriched genes derived from large-scale transcriptomics studies. 

Simulation studies show that BioQC is both fast and sensitive in detecting tissue heterogeneity. 

## Usage

The usage is documented on our [website](https://accio.github.io/BioQC) or in the following two
vignettes: 

* [Basic Usage](vignettes/bioqc.md)
* [Gene Set Analysis with BioQC](vignettes/bioqc-signedGenesets.md)


## Case studies

There are case-studies available on our [website](https://accio.github.io/BioQC) or in the [BioQC-example](https://github.com/Accio/BioQC-example) github repository. 

## Citation

Please cite the following publication if you find BioQC useful:

Zhang, Jitao David, Klas Hatje, Gregor Sturm, Clemens Broger, Martin Ebeling, Martine Burtin, Fabiola Terzi, Silvia Ines Pomposiello, and Laura Badi. “Detect Tissue Heterogeneity in Gene Expression Data with BioQC.” BMC Genomics 18 (2017): 277. doi:10.1186/s12864-017-3661-2.

The open-access text of the publication can be found [on the website of BMC Genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3661-2).
