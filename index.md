---
title: Detect tissue heterogeneity in gene expression data with BioQC
sidebar: bioqc_sidebar
permalink: index.html
---

BioQC is a R/[Bioconductor](https://bioconductor.org/packages/release/bioc/html/BioQC.html) package to detect tissue heterogeneity in gene expression data. To this end BioQC implements a computationally efficient Wilcoxon-Mann-Whitney test and provides more than 150 signatures of tissue-enriched genes derived from large-scale transcriptomics studies. 

Simulation studies show that BioQC is both fast and sensitive in detecting tissue heterogeneity. 

BioQC is part of the following publication: 

> Jitao David Zhang, Klas Hatje, Clemens Broger, Martin Ebeling, Martine Burtin, Fabiola Terzi, Silvia Ines Pomposiello, Gregor Sturm and Laura Badi "Detect tissue heterogeneity in gene expression data with BioQC", *BMC Genomics*: [*submitted*]. 

## Usage
In this section, we demonstrate how to apply BioQC to a gene expression dataset. Moreover, we 
show that it is also feasible to perform a gene set enrichment analysis with BioQC. 
* [Basic Usage](bioqc.html)
* [Gene Set Analysis with BioQC](bioqc-signedGenesets.html)

## The Algorithm
We explain why we chose the Wilcoxon-Mann-Whitney test as underlying statistical method and explain the algorithmic improvments we implemented to make it feasible to apply BioQC on large-scale studies. We show, that our optimized version of the Wilcoxon-Mann-Whitney test outperforms the native R implementation by more than 1000-times.  
* [WMW test comparison](bioqc-wmw-test-performance.html)
* [Algorithmic improvements](bioqc-efficiency.html)

## Case studies
We present two case-studies, the first is a practical example where BioQC successfully identified contamination of pancreas in a liver gene expression study. In the second, we benchmark sensitivity and specificity of BioQC with both simulated and real data. 
* [Kidney Expression Example](bioqc-kidney.html)
* [Benchmark: Sensitivity and Specificity](bioqc-simulation.html)

{% include links.html %}
