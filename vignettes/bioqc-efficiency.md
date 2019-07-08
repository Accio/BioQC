---
title: "BioQC Algorithm: Speeding up the Wilcoxon-Mann-Whitney Test"
author: "Gregor Sturm and Jitao David Zhang"
package: BioQC
date: "2019-07-08"
vignette: >
  %\VignetteIndexEntry{BioQC Alogrithm: Speeding up the Wilcoxon-Mann-Whitney Test}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
output: 
  rmarkdown::html_vignette:
    self_contained: no
  md_document:
    variant: markdown_phpextra
    preserve_yaml: TRUE
---

Supplementary Information for "Detect issue heterogenity in gene expression data with [*BioQC*](https://github.com/Accio/BioQC)" ([Jitao David Zhang](mailto:jitao_david.zhang@roche.com), Klas Hatje, Gregor Sturm, Clemens Broger, Martin Ebeling, Martine Burtin, Fabiola Terzi, Silvia Ines Pomposiello and [Laura Badi](laura.badi@roche.com))




In this vignette, we explain the underlaying algorithmic details of our implementation of the Wilcoxon-Mann-Whitney test. The source code used to produce this document can be found in the github repository [BioQC](https://github.com/Accio/BioQC/vignettes).

*BioQC* is a R/Bioconductor package to detect tissue heterogeneity from high-throughput gene expression profiling data. It implements an efficient Wilcoxon-Mann-Whitney test, and offers tissue-specific gene signatures that are ready to use 'out of the box'.













