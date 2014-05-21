\name{wmwTest-methods}
\docType{methods}
\alias{wmwTest}
\alias{wmwTest-methods}
\alias{wmwTest,ExpressionSet,list,character,logical-method}
\alias{wmwTest,ExpressionSet,logical,character,logical-method}
\alias{wmwTest,ExpressionSet,numeric,character,logical-method}
\alias{wmwTest,matrix,logical,character,logical-method}
\alias{wmwTest,matrix,numeric,character,logical-method}
\alias{wmwTest,numeric,logical,character,logical-method}
\alias{wmwTest,numeric,numeric,character,logical-method}
\alias{wmwTest,ANY,ANY,missing,ANY-method}
\alias{wmwTest,ANY,ANY,ANY,missing-method}
\alias{wmwTest,ANY,ANY,missing,missing-method}
\title{Wilcoxon-Mann-Whitney test for gene signatures}
\description{
  The Wilcoxon-Mann-Whitney Test (WMW-Test) is a non-parameteric statistic tool to examine whether one of the two samples of independent observations tends to have larger values than the other. It can be used in cases where higher (and/or lower) rankings of a group of genes are of interest. We illustrate the use of the function with the example of using tissue-specific marker genes to predict source organ of microarray samples.

  The basic form of \emph{wmwTest} function takes a numeric vector (\code{object}) and a  boolean vector (\code{sub}) indicating which genes belong to the markers (\code{TRUE}) and which not (\code{FALSE}). In the context of our example, the numeric vector contains expression values of all genes in a sample, and the \code{TRUE} values in \code{sub} specify tissue-specific marker genes. Then the rank of these tissue-specific genes are compared to the rest of genes using the Wilcoxon-Mann-Whitney test. If they are generally higher expressed than other genes, there is higher possibility that the tested sample may come from the same tissue where the marker genes were derived. Running the test over a set of tissue gene markers, and comparing statistics or p-values, then provides a naive way to predict the most probable tissue(s) from which the sample was taken.

  In case only the over-presentation of genes is of interest (as in the
  tissue-heterogeneity detection example), one-sided alternative
  hypothesis \code{greater} should be used. The
  \code{two.sided} hypothesis or the other one-sided hypothesis
  (\code{less}) can be used depending on the purpose of test.
  
  The function returns the \code{p value} of the MWM-test if the option \code{statistic} is set to \code{FALSE}, otherwise the statistic \code{W} is returned.

  Other forms of the methods are provided as wrappers for more
  complex data structures.
}
\section{Methods}{
\describe{
  \item{\code{signature(object = "numeric", sub = "logical", alternative = "character", statistic = "logical")}}{
    Basic form of the function. See descriptions above.
    
    The \code{alternative} and \code{statistic} options can be missing
    (although not recommended): \code{alternative} has the default value \code{two.sided}, and
    \code{statistic} the default value \code{FALSE}.
  }
  
  \item{\code{signature(object = "numeric", sub = "numeric", alternative = "character", statistic = "logical")}}{
    \code{sub} contains numeric indices (probably integers) indicating
    which positions are specified as markers.
  }
  
  \item{\code{signature(object = "matrix", sub = "logical", alternative = "character", statistic = "logical")}}{
    The WMW-Test is run with the same set of logical vectors on all
    columns of the matrix.
  }
  \item{\code{signature(object = "matrix", sub = "numeric", alternative = "character", statistic = "logical")}}{
    The WMW-Test is run with the same set of numerical vectors on all
    columns of the matrix. The numerical vectors specify which rows are
    to be compared with rest of the rows.
  }
  
  \item{\code{signature(object = "ExpressionSet", sub = "logical", alternative = "character", statistic = "logical")}}{
    The WMW-Test is run on each sample of the \code{ExpressionSet}
    object with the same set of logical vectors. \code{TRUE} value of
    the vector specify which features to be compared to rest features.
  }
  \item{\code{signature(object = "ExpressionSet", sub = "numeric", alternative = "character", statistic = "logical")}}{
    The WMW-Test is run on each sample of the \code{ExpressionSet}
    object with the same set of numeric vectors. The numeric vector
    specify which features to be compared to the rest of features.
  }

  \item{\code{signature(object = "ExpressionSet", sub = "list", alternative = "character", statistic = "logical")}}{
    The \code{sub} list can either be a list of numeric or logical
    vectors, and in such cases the function checks each signature by
    calling \code{lapply}.

    More interestingly, the list can be returned by the \code{read_gmt_list}
    function in the \code{ribiiosIO} package. In such \code{GMT} lists, each item of the list must at least contain a string named \code{name}, and a
    character vector named \code{genes}. In this approach signatures
    stored in GMT files can be easily checked for over- or
    under-expression by using the WMW-Test. To this purpose a column named
    \code{GeneSymbol} must be present in the \code{fData} of
    \code{ExpressionSet}: the values of these column will be matched to
    the gene signatures in the \code{name} vector of the \code{sub}
    list.

    Samples are examined with all signature genes in the \code{sub} list, which most probably come from a \code{GMT}
    file. The function returns a matrix with tissue-markers as rows and
    samples as columns. Each sample has a p-value (or W-statistic)
    associated with each tissue signature gene set.}
}
}
\examples{
## basic forms
rd <- rnorm(1000)
rl <- sample(c(TRUE, FALSE), 1000, replace=TRUE)
wmwTest(rd, rl, alternative="two.sided", statistic=FALSE)
wmwTest(rd, which(rl), alternative="two.sided", statistic=FALSE)
rd1 <- rd + ifelse(rl, 0.5, 0)
wmwTest(rd1, rl, alternative="greater", statistic=FALSE)
wmwTest(rd1, rl, alternative="greater", statistic=TRUE)
rd2 <- rd - ifelse(rl, 0.2, 0)
wmwTest(rd2, rl, alternative="greater", statistic=FALSE)
wmwTest(rd2, rl, alternative="two.sided", statistic=FALSE)
wmwTest(rd2, rl, alternative="less", statistic=FALSE)

## matrix forms
rmat <- matrix(c(rd, rd1, rd2), ncol=3, byrow=FALSE)
wmwTest(rmat, rl, alternative="two.sided", statistic=FALSE)
wmwTest(rmat, rl, alternative="greater", statistic=FALSE)

## using ExpressionSet
if(require(Biobase)) {
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

  ## TODO: using GMT lists
}
}
\keyword{methods}