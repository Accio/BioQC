---
title: "Comparing the Wilcoxon-Mann-Whitney to alternative statistical tests"
permalink: bioqc-wmw-test-performance.html
output: 
  md_document:
    variant: markdown_phpextra
    preserve_yaml: TRUE
---

In this document, we show that the Wilcoxon-Mann-Whitney test is
comparable or superior to alternative methods.

Two alternative methods could be compared with the Wilcoxon-Mann-Whitney
(WMW) test proposed by BioQC: the Kolmogorov-Smirnov (KS) test, and the
Student’s t-test, or more particularly, the Welch’s test which does not
assume equal sample number or equal variance, which is appropriate in
the setting of gene expression studies.

1.  It is documented in statistics literature that the WMW test offers a
    higher power than the Kolmogorov-Smirnov test[^1],[^2].
2.  Compared with parameterized test methods such as the t-test, the WMW
    test is (a) resistance to monotone transformation, (b) suffers less
    from outliers, and (c) provides higher efficiency when many genes
    are profiled and the distribution of gene expression deviates from
    the normal distribution, which are important criteria in genome-wide
    expression data.

Based on these considerations, BioQC implements a computationally
efficient version of the WMW test. In order not to confuse end-users, no
alternative methods are implemented.

Nevertheless, in order to demonstrate the power of WMW test in
comparison with the KS-test or the t-test, we performed the sensitivity
benchmark described in the [simulation studies](bioqc-simulation.html),
for the two alternative tests respectively.

<img src="pages/bioqc/bioqc-wmw-test-performance_files/figure-markdown_phpextra/sensitivity_benchmark_fig-1.svg" alt="**Figure 1:** Sensitivity benchmark. Expression levels of genes in the ovary signature are dedicately sampled randomly from normal distributions with different mean values. The lines show the enrichment score for the Wilcoxon-Mann-Whitney test, the t-test and the Kolmogorov-Smirnov test respectively. In the right panel, outliers were added by adding a random value to 1% of the simulated genes. " style="display:block; margin: auto" />
<p markdown="1" class="caption">
**Figure 1:** Sensitivity benchmark. Expression levels of genes in the
ovary signature are dedicately sampled randomly from normal
distributions with different mean values. The lines show the enrichment
score for the Wilcoxon-Mann-Whitney test, the t-test and the
Kolmogorov-Smirnov test respectively. In the right panel, outliers were
added by adding a random value to 1% of the simulated genes.
</p>

As expected, the results suggest, that both the KS-test and the WMW-test
are robust to noise, while the performance of the t-test drops
significantly on noisy data. Additionally, the WMW-test appears to be
superior to the KS-test for low expression differences.

Computational Performance {#computational-performance}
-------------------------

Since the KS-test is so slow, we did not replicate the sensitivity
benchmark from the [simulation studies](bioqc-simulation.html) using the
enrichment score rank. While it takes BioQC about 6 seconds on a single
thread to test all 155 signatures, it already takes the KS-test about 8
seconds to test a single signature.

    ##       test replications elapsed relative
    ## 2  runKS()           10  75.200     1.16
    ## 1 runWMW()           10  64.816     1.00

R Session Info {#r-session-info}
--------------

~~~~ r
sessionInfo()
~~~~

    ## R version 3.3.1 (2016-06-21)
    ## Platform: i686-pc-linux-gnu (32-bit)
    ## Running under: Linux Mint 18
    ## 
    ## locale:
    ##  [1] LC_CTYPE=de_CH.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=de_CH.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_CH.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
    ## [8] datasets  base     
    ## 
    ## other attached packages:
    ##  [1] plyr_1.8.4           rbenchmark_1.0.0     ggplot2_2.2.0       
    ##  [4] reshape2_1.4.2       gplots_3.0.1         xtable_1.8-2        
    ##  [7] GEOquery_2.36.0      gridExtra_2.2.1      latticeExtra_0.6-28 
    ## [10] RColorBrewer_1.1-2   lattice_0.20-34      hgu133plus2.db_3.2.2
    ## [13] org.Hs.eg.db_3.2.3   RSQLite_1.0.0        DBI_0.5-1           
    ## [16] AnnotationDbi_1.32.3 IRanges_2.4.8        S4Vectors_0.8.11    
    ## [19] BioQC_1.02.1         Biobase_2.30.0       BiocGenerics_0.16.1 
    ## [22] Rcpp_0.12.8          testthat_1.0.2       knitr_1.15          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gdata_2.17.0       bitops_1.0-6       stringr_1.1.0     
    ##  [4] gtools_3.5.0       labeling_0.3       highr_0.6         
    ##  [7] tools_3.3.1        RCurl_1.95-4.8     rmarkdown_1.1     
    ## [10] R6_2.2.0           scales_0.4.1       assertthat_0.1    
    ## [13] digest_0.6.10      gtable_0.2.0       evaluate_0.10     
    ## [16] KernSmooth_2.23-15 stringi_1.1.2      caTools_1.17.1    
    ## [19] htmltools_0.3.5    munsell_0.4.3      grid_3.3.1        
    ## [22] XML_3.98-1.5       colorspace_1.3-1   magrittr_1.5      
    ## [25] lazyeval_0.2.0     yaml_2.1.13        crayon_1.3.2      
    ## [28] tibble_1.2

References {#references}
----------

[^1]: Irizarry, Rafael A., et al. "Gene set enrichment analysis made
    simple."Statistical methods in medical research 18.6 (2009):
    565-575.

[^2]: Filion, Guillaume J. "The signed Kolmogorov-Smirnov test: why it
    should not be used." GigaScience 4.1 (2015): 1.
