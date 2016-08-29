<!--
%\VignetteEngine{knitr::docco_linear}
%\VignetteIndexEntry{Using BioQC with signed genesets}
-->

Using BioQC with signed genesets
================================


Introduction
------------

Besides being used as a QC-tool for gene expression data, **BioQC** can be used for general-purpose gene-set enrichment analysis. Its core function, _wmwTest_, can handle not only "normal" gene sets, but also signed genesets where two sets of genes represent signatures that are positively and negatively regulated, respectively. 

This vignette describes an example of such applications.

GMT read-in
----------
We first open a toy GMT file containing signed genesets with _readSignedGmt_. The function reads the file into a _SignedGenesets_ object.


```r
library(BioQC)
gmtFile <- system.file("extdata/test.gmt", package="BioQC")
## print the file context
cat(readLines(gmtFile), sep="\n")
```

```
## GS_A	TestGeneSet_A	AKT1	AKT2	AKT3
## GS_B	TestGeneSet_B	MAPK1	MAPK3	MAPK8
## GS_C_UP	TestGeneSet_C	ERBB2	ERBB3
## GS_C_DN	TestGeneSet_C	EGFR	ERBB4
## GS_D_UP	TestGeneSet_D	GATA2	GATA4
## GS_D_DN	TestGeneSet_D	GATA1	GATA3
## GS_E_DN	TestGeneSet_E	TSC1	TSC2
```

```r
## read the file into SignedGenesets
genesets <- readSignedGmt(gmtFile)
print(genesets)
```

```
## A list of 5 signed gene-sets
##   GS_A (no genes)
##   GS_B (no genes)
##   GS_C[n=4]
##     positive[n=2]:ERBB2,ERBB3
##     negative[n=2]:EGFR,ERBB4
##   GS_D[n=4]
##     positive[n=2]:GATA2,GATA4
##     negative[n=2]:GATA1,GATA3
##   GS_E[n=2]
##     positive[n=0]:
##     negative[n=2]:TSC1,TSC2
```

Note that though the GMT file contains 6 non-empty lines, only five signed genesets are returned by _readSignedGmt_ because a pair of negative and negative geneset of GS_C is available.

Note that in its current form, no genes are reported for GS_A and GS_B, because the names of both sets do not match the patterns. User can decide how to handle such genesets that match patterns of neither positive nor negative sets (GS_A and GS_B in this case). By default, such genesets are ignored; however, user can decide to treat them as either positive or negative genesets. For the purpose of this tutorial, we treat these genesets as positive.


```r
genesets <- readSignedGmt(gmtFile, nomatch="pos")
print(genesets)
```

```
## A list of 5 signed gene-sets
##   GS_A[n=3]
##     positive[n=3]:AKT1,AKT2,AKT3
##     negative[n=0]:
##   GS_B[n=3]
##     positive[n=3]:MAPK1,MAPK3,MAPK8
##     negative[n=0]:
##   GS_C[n=4]
##     positive[n=2]:ERBB2,ERBB3
##     negative[n=2]:EGFR,ERBB4
##   GS_D[n=4]
##     positive[n=2]:GATA2,GATA4
##     negative[n=2]:GATA1,GATA3
##   GS_E[n=2]
##     positive[n=0]:
##     negative[n=2]:TSC1,TSC2
```

Expression object construction
------------------------------
Next we construct an ExpressionSet object containing 100 genes, with some of the genes in the GMT file present. 


```r
set.seed(1887)
testN <- 100L
testSampleCount <- 3L
testGenes <- c("AKT1", "AKT2", "ERBB2", "ERBB3", "EGFR","TSC1", "TSC2", "GATA2", "GATA4", "GATA1", "GATA3")
testRows <- c(testGenes, paste("Gene", (length(testGenes)+1):testN, sep=""))
testMatrix <- matrix(rnorm(testN*testSampleCount, sd=0.1), nrow=testN, dimnames=list(testRows, NULL))
testMatrix[1:2,] <- testMatrix[1:2,]+10
testMatrix[6:7,] <- testMatrix[6:7,]-10
testMatrix[3:4,] <- testMatrix[3:4,]-5
testMatrix[5,] <- testMatrix[5,]+5
testEset <- new("ExpressionSet", exprs=testMatrix)
```
The expression of somes genes are delibrately changed so that some genesets are expected to be significantly enriched.

* AKT1 and AKT2 are higher expressed than the background genes ('background' hereafter). 
* TSC1 and TSC2 are lower expressed than the background.
* ERBB2 and ERBB3 are moderately lower expressed than the background.
* EGFR is moderately higher expressed than the background.

Considering the geneset composition printed above, we can expect that
* GS_A shall be significantly up-regulated; 
* GS_B shall not be significant, because its genes are missing;
* GS_C shall be significantly down-regulated, with positive targets down-regulated and negative targets up-regulated;
* GS_D shall not always be significant, because its genes are comparable to the background;
* GS_E shall be significantly up-regulated, because its negative targets are down-regulated.


Match genes
-----------

With both _SignedGenesets_ and _ExpressionSet_ objects at hand, we can query the indices of genes of genesets in the _ExpressionSet_ object.


```r
testIndex <- matchGenes(genesets, testEset, col=NULL)
print(testIndex)
```

```
## A list of 5 signed indices with offset=1
## Options: NA removed: TRUE; duplicates removed: TRUE
##   NONAME[n=2]
##     positive[n=2]:1,2
##     negative[n=0]:
##   NONAME (no genes)
##   NONAME[n=3]
##     positive[n=2]:3,4
##     negative[n=1]:5
##   NONAME[n=4]
##     positive[n=2]:8,9
##     negative[n=2]:10,11
##   NONAME[n=2]
##     positive[n=0]:
##     negative[n=2]:6,7
```

If the _ExpressionSet_ object has GeneSymbols in its featureData other than the feature names, user can specify the parameter _col_ in _matchGenes_ to make the match working. _matchGene_ also works with characters and matrix as input, please check out the help pages for such uses.

Perform the analysis
--------------------


```r
wmwResult.greater <- wmwTest(testEset, testIndex, valType="p.greater")
print(wmwResult.greater)
```

```
##                1           2           3
## [1,] 0.008185751 0.008185751 0.008185751
## [2,] 1.000000000 1.000000000 1.000000000
## [3,] 0.997942914 0.997942914 0.997942914
## [4,] 0.924153922 0.114118870 0.211756278
## [5,] 0.008185751 0.008185751 0.008185751
```

As expected, GS\_A and GS\_E are significant when the alternative hypothesis is _greater_.


```r
wmwResult.less <- wmwTest(testEset, testIndex, valType="p.less")
print(wmwResult.less)
```

```
##                1           2           3
## [1,] 0.992348913 0.992348913 0.992348913
## [2,] 1.000000000 1.000000000 1.000000000
## [3,] 0.002192387 0.002192387 0.002192387
## [4,] 0.078389208 0.889240837 0.793302033
## [5,] 0.992348913 0.992348913 0.992348913
```
As expected, GS\_C is significant when the alternative hypothesis is _less_.


```r
wmwResult.two.sided <- wmwTest(testEset, testIndex, valType="p.two.sided")
print(wmwResult.two.sided)
```

```
##                1           2           3
## [1,] 0.016371502 0.016371502 0.016371502
## [2,] 1.000000000 1.000000000 1.000000000
## [3,] 0.004384774 0.004384774 0.004384774
## [4,] 0.156778415 0.228237741 0.423512555
## [5,] 0.016371502 0.016371502 0.016371502
```

As expected, GS\_A, GS\_C, and GS\_E are significant when the alternative hypothesis is _two.sided_. 

A simple way to integrate the results of both alternative hypotheses into one numerical value is the _Q_ value, which is defined as absolute log10 transformation of the smaller p value of the two hypotheses, times the sign defined by _greater_=1 and _less_=-1.


```r
wmwResult.Q <- wmwTest(testEset, testIndex, valType="Q")
print(wmwResult.Q)
```

```
##               1          2          3
## [1,]  1.7859115  1.7859115  1.7859115
## [2,]  0.0000000  0.0000000  0.0000000
## [3,] -2.3580528 -2.3580528 -2.3580528
## [4,] -0.8047137  0.6416125  0.3731337
## [5,]  1.7859115  1.7859115  1.7859115
```

As expected, genesets GS\_A, GS\_C, and GS\_E have large absolute values, suggesting they are potentially enriched; while GS\_A and GS\_E are positively enriched, GS\_C is negatively enriched.

Conclusions
-----------
In summary, _BioQC_ can be used to perform gene-set enrichment analysis with signed genesets. Its speed is comparable to that of unsigned genesets.

Acknowledgement
----------------
I would like to thank the authors of the gCMAP package with the idea of integrating signed genesets into the framework of Wilcoxon-Mann-Whitney test.

R session info
--------------

```r
sessionInfo()
```

```
## R version 3.3.0 (2016-05-03)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux Server release 6.3 (Santiago)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] knitr_1.13          BioQC_1.02.0        Biobase_2.32.0     
## [4] BiocGenerics_0.18.0 Rcpp_0.12.5        
## 
## loaded via a namespace (and not attached):
## [1] compiler_3.3.0 magrittr_1.5   formatR_1.4    markdown_0.7.7
## [5] tools_3.3.0    stringi_1.1.1  stringr_1.0.0  evaluate_0.9
```
