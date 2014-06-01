\name{wmwTest}
\alias{wmwTest}
\title{
  Wilcoxon-Mann-Whitney rank sum test for high-throughput expression
  profiling data
}
\description{
  We have implemented an highly efficient Wilcoxon-Mann-Whitney rank sum
  test for high-throughput expression profiling data. For datasets with
  more than 100 features (genes), the function performs almost identical
  to its R implementations (\code{wilcox.test} in \code{stats}, or
  \code{rankSumTestWithCorrelation} in \code{limma}) can be more than
  1000 times faster.
}
\usage{
wmwTest(x, ind.list, alternative = c("greater", "less", "two.sided", "U"), simplify = TRUE)
}
\arguments{
  \item{x}{A numeric matrix. All other data types (e.g. numeric vectors
    or \code{ExpressionSet} objects) are coerced into matrix.}
  \item{ind.list}{A list of integer indices (starting from 1) indicating
  signature genes. Can be of length zero. Other data types (e.g. a list
  of numeric or logical vectors, or a numeric or logical vector) are
  coerced into such a list. See \code{details} below for a special case
  using GMT files.}
  \item{alternative}{The value type to be returned, allowed values
    include \code{greater} and \code{less} (one-sided tests),
    \code{two.sided}, and \code{U} statistic.}
  \item{simplify}{Logical. If not, the returning value is in matrix
    format; if set to \code{TRUE}, the results are simplified into
    vectors when possible (default).}
}
\details{
  The basic application of the function is to test the enrichment of
  gene sets in expression profiling data or differentially expressed
  data (the matrix with feature/gene in rows and samples in columns).

  A special case is when \code{x} is an \code{eSet} object
  (e.g. \code{ExpressionSet}), and \code{ind.list} is a list returned
  from \code{readGmt} function. In this case, the only requirement is
  that one column named \code{GeneSymbol} in the \code{featureData}
  contain gene symbols used in the GMT file. See the example below.
}
\value{
  A numeric matrix or vector containing the statistic.
}
\references{
Barry, W.T., Nobel, A.B., and Wright, F.A. (2008). A statistical framework for testing functional categories in microarray data. _Annals of Applied Statistics_ 2, 286-315.

Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test
accounting for inter-gene correlation. _Nucleic Acids Research_ 40(17):e133

Zar, JH (1999). _Biostatistical Analysis 4th Edition_. Prentice-Hall International, Upper Saddle River, New Jersey.
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\note{
  The function has been optimized for expression profiling data. It
  avoids repetitive ranking of data as done by native R implementations
  and uses efficient C code to increase the performance and control
  memory use. Simulation studies using expression profiles of 22000
  genes in 2000 samples and 200 gene sets suggested that the C
  implementation can be >1000 times faster than the R
  implementation. And it is possible to further accelerate by
  parallel calling the function with \code{mclapply} in the \code{multicore} package.
}

\seealso{
  \code{wilcox.test} in the \code{stats} package, and \code{rankSumTestWithCorrelation} in
  the \code{limma} package.
}
\examples{
## R-native data structures
set.seed(1887)
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

wmwTest(rmat, which(rl), alternative="two.sided")
wmwTest(rmat, which(rl), alternative="greater")

## using ExpressionSet
data(sample.ExpressionSet)
testSet <- sample.ExpressionSet
fData(testSet)$GeneSymbol <- paste("GENE_",1:nrow(testSet), sep="")
mySig1 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.25, 0.75), replace=TRUE)
wmwTest(testSet, which(mySig1), alternative="greater")

## using integer
exprs(testSet)[,1L] <- exprs(testSet)[,1L] + ifelse(mySig1, 50, 0)
wmwTest(testSet, which(mySig1), alternative="greater")

## using lists
mySig2 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.6, 0.4), replace=TRUE)
wmwTest(testSet, list(first=mySig1, second=mySig2))

## using GMT file
gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt_list <- readGmt(gmt_file)

gss <- sample(unlist(sapply(gmt_list, function(x) x$genes)), 1000)
eset<-new("ExpressionSet",
          exprs=matrix(rnorm(10000), nrow=1000L),
          phenoData=new("AnnotatedDataFrame", data.frame(Sample=LETTERS[1:10])),
          featureData=new("AnnotatedDataFrame",data.frame(GeneSymbol=gss)))
esetWmwRes <- wmwTest(eset ,gmt_list, alternative="greater")
summary(esetWmwRes)
}