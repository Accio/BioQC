\name{wmwTest-methods}
\docType{methods}
\alias{wmwTest}
\alias{wmwTest-methods}
\alias{wmwTest,ExpressionSet,list,character-method}
\alias{wmwTest,ExpressionSet,gmtlist,character-method}
\alias{wmwTest,matrix,list,character-method}
\alias{wmwTest,numeric,numeric,character-method}
\alias{wmwTest,numeric,logical,character-method}
\alias{wmwTest,matrix,numeric,character-method}
\alias{wmwTest,matrix,logical,character-method}
\alias{wmwTest,ANY,ANY,missing-method}
\title{Wilcoxon-Mann-Whitney test for gene signatures}
\description{
  The Wilcoxon-Mann-Whitney rank sum Test (WMW-Test) is a
  non-parameteric statistic tool to examine whether one of the two
  samples of independent observations tends to have larger values than
  the other. It can be used in cases where higher (and/or lower)
  rankings of a group of genes are of interest. We apply the this test
  to detect tissue heterogeneity in gene expression data using
  tissue-specific marker genes. The statistical assumptions include (1)
  gene expression levels are independent from each other, (2) expression
  profiling readouts is rank-invariant with regard to true expression
  data, and (3) tissue contanmination (either due to technical or
  biological reasons) lead to significant increase of tissue-specific
  gene's expresion. Some of these hypothese are very strong, however we
  found in reality the statistical test, despite of its simplicity,
  works wells and often detects potentia issues in expression profiling
  data.

  The basic form of \emph{wmwTest} function takes a numeric matrix and a
  list of integers (starting from 1 following R's convention) indicating
  which genes belong to the gene set. All other data structures are cast
  to this data structure, since the algorithm has been optimized for use with
  expression data.

  In case only the over-presentation of genes is of interest (as in the
  tissue-heterogeneity detection example), one-sided alternative
  hypothesis \code{greater} should be used (which is also the default
  setting if \code{alternative} is missing). The
  \code{two.sided} hypothesis or the other one-sided hypothesis
  (\code{less}) can be used depending on the purpose of test.

  Other forms of the methods are provided as wrappers for more
  complex data structures.
}
\section{Methods}{
\describe{
  \item{\code{signature(exprs = "matrix", index = "numeric", alternative = "character")}}{
    Basic form of the function. Rows of the matrix contain genes, and
    columns contain samples.

    Though \code{index} should be integers starting from 1, numeric values are tolerated
    and automatically casted.
  }
  
  \item{\code{signature(exprs = "matrix", index = "logical", alternative = "character")}}{
    \code{index} contains logical vectors of the same length as \code{nrow(exprs)} indicating
    which positions are specified as markers (\code{TRUE}).
  }
  
  \item{\code{signature(exprs = "numeric", index = "logical", alternative = "character")}}{
    The WMW test est is run with one numeric vector.
  }
  \item{\code{signature(exprs = "numeric", index = "numeric", alternative = "character")}}{
    The WMW test est is run with one numeric vector.
  }
  
  \item{\code{signature(exprs = "matrix", index = "list", alternative = "character")}}{
    The WMW test is run on a series of samples with a list of either
    integer or logical values indicating which rows are marker genes.
  }
  \item{\code{signature(exprs = "eSet", index = "list", alternative = "character")}}{
    The WMW test is run on an \code{eSet} object
    (e.g. \code{ExpressionSet}) with a list of either
    integer or logical values indicating which rows are marker genes.
  }
  \item{\code{signature(exprs = "eSet", index = "gmtlist", alternative = "character")}}{
    \code{index} is be returned by the \code{read_gmt_list}
    function, which reads in (even large) GMT file efficiently.

    In such \code{GMT} lists, each item of the list must at least contain a string named \code{name}, and a
    character vector named \code{genes}. In this approach signatures
    stored in GMT files can be easily checked for over- or
    under-expression by using the WMW-Test.

    A column named \code{GeneSymbol} must be present in the \code{fData} of
    \code{eSet}: the values of this column will be matched to
    the gene signatures in the \code{name} vector of the \code{index}
    list.}
}
}
\examples{
## basic forms
rd <- rnorm(1000)
rl <- sample(c(TRUE, FALSE), 1000, replace=TRUE)
wmwTest(rd, rl, alternative="two.sided")
wmwTest(rd, which(rl), alternative="two.sided")
rd1 <- rd + ifelse(rl, 0.5, 0)
wmwTest(rd1, rl, alternative="greater")
wmwTest(rd1, rl, alternative="U")
rd2 <- rd - ifelse(rl, 0.2, 0)
wmwTest(rd2, rl, alternative="greater")
wmwTest(rd2, rl, alternative="two.sided")
wmwTest(rd2, rl, alternative="less")

## matrix forms
rmat <- matrix(c(rd, rd1, rd2), ncol=3, byrow=FALSE)
wmwTest(rmat, rl, alternative="two.sided")
wmwTest(rmat, rl, alternative="greater")

## using ExpressionSet
data(sample.ExpressionSet)
testSet <- sample.ExpressionSet
fData(testSet)$GeneSymbol <- paste("GENE_",1:nrow(testSet), sep="")
mySig1 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.25, 0.75), replace=TRUE)
wmwTest(testSet, mySig1, alternative="greater",statistic=FALSE)

## using integer
exprs(testSet)[,1L] <- exprs(testSet)[,1L] + ifelse(mySig1, 50, 0)
wmwTest(testSet, which(mySig1), alternative="greater",statistic=FALSE)

## using lists
mySig2 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.6, 0.4), replace=TRUE)
wmwTest(testSet, list(first=mySig1, second=mySig2))

## using GMT file
}
\keyword{methods}