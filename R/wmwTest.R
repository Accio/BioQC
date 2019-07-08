#' @useDynLib BioQC, .registration=TRUE, .fixes="C_"
#' @importFrom Rcpp sourceCpp

## FOR REFERENCE rankSumTestWithCorrelation function from the limma package (version 3.18.13)

## authors: Gordon Smyth and Di Wu, following Zar, JD (1999) Biostatistical Analysis 4th Edition
## used under GPL(>=2) license. The function has been slihtly modified to allow reporting results

#' @importFrom stats pt
rankSumTestWithCorrelation <- function (index, statistics, correlation = 0, df = Inf) {
    n <- length(statistics)
    r <- rank(statistics)
    r1 <- r[index]
    n1 <- length(r1)
    n2 <- n - n1
    U <- n1 * n2 + n1 * (n1 + 1)/2 - sum(r1)
    mu <- n1 * n2/2
    if (correlation == 0 || n1 == 1) {
        sigma2 <- n1 * n2 * (n + 1)/12
    }
    else {
        sigma2 <- asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 - 
            1) + asin(correlation/2) * n1 * (n1 - 1) * n2 * (n2 - 
            1) + asin((correlation + 1)/2) * n1 * (n1 - 1) * 
            n2
        sigma2 <- sigma2/2/pi
    }
    TIES <- (length(r) != length(unique(r)))
    if (TIES) {
        NTIES <- table(r)
        prod <- sum(NTIES * (NTIES + 1) * (NTIES - 1))
        denom <- (n * (n + 1) * (n - 1))
        adjustment <- prod/denom
        sigma2 <- sigma2 * (1 - adjustment)
    }
    zlowertail <- (U + 0.5 - mu)/sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu)/sqrt(sigma2)
    less <-pt(zuppertail, df = df, lower.tail = FALSE)
    greater <- pt(zlowertail, df = df)
    res <- c(U=U,
             mu=mu,
             n1=n1,
             n2=n2,
             sigma2=sigma2,
             r1sum=sum(r1),
             zlt=zlowertail,
             zut=zuppertail,
             less = less,
             greater = greater)
    return(res)
}


#' Wilcoxon-Mann-Whitney test in R
#'
#' @param x A numerical vector
#' @param sub A logical vector or integer vector to subset \code{x}. Numbers in \code{sub} are compared with numbers out of \code{sub}
#' @param valType Type of retured-value. Supported values: p.greater, p.less, p.two.sided, and W statistic (note it is different from the U statistic)
#' 
#' @importFrom stats wilcox.test
#' 
#' @examples
#' testNums <- 1:10
#' testSub <- rep_len(c(TRUE, FALSE), length.out=length(testNums))
#' wmwTestInR(testNums, testSub)
#' wmwTestInR(testNums, testSub, valType="p.two.sided")
#' wmwTestInR(testNums, testSub, valType="p.less")
#' wmwTestInR(testNums, testSub, valType="W")
#' @export
wmwTestInR <- function(x, sub, valType=c("p.greater", "p.less", "p.two.sided",  "W")) {
    if(is.numeric(sub)) {
        tmp <- rep(FALSE, length(x))
        tmp[sub] <- TRUE
        sub <- tmp
    }
    valType <- match.arg(valType)
    if(!is.logical(sub))
        stop("sub must be either numeric indices or logical")
    isStat <- valType=="W"
    if(!any(sub)) return(ifelse(isStat, 0, 1))
    if(!isStat) {
        alternative <- substr(valType, 3, nchar(valType))
        stopifnot(alternative %in% c("greater", "less", "two.sided"))
    } else {
        alternative <- "two.sided"
    }
    wt <- wilcox.test(x[sub],x[!sub],
                      alternative=alternative, exact=FALSE)
    return(ifelse(isStat, wt$statistic, wt$p.value))
}

## type2int and formatMatrixInd are helper functions for wmwTest and wmwSignedTest
TYPE_CODES <- c("p.greater"=0L, "p.less"=1L,
                "p.two.sided"=2L, "U"=3L,
                "abs.log10p.greater"=4L,
                "log10p.less"=5L,
                "abs.log10p.two.sided"=6L,
                "Q"=7L)
#'prints the options of valTypes of wmwTest
valTypes <- function() names(TYPE_CODES)

type2int <- function(type) {
    if(!type %in% names(TYPE_CODES))
        stop("Should not happen! This is the wrong code")
    return(TYPE_CODES[type])
}

##----------------------------------------##
## wmwTest
##----------------------------------------##
#' @importFrom methods is
wmwTest.default <- function(matrix,
                            indexList,
                            valType=c("p.greater", "p.less", "p.two.sided", "U",
                                "abs.log10p.greater","log10p.less","abs.log10p.two.sided",
                                "Q"),
                            simplify=TRUE) {
    if(!is.matrix(matrix) || !is(indexList, "IndexList"))
        stop("'matrix' and 'indexList' must be matrix and an IndexList object, respectively")
    if(missing(simplify))  simplify <- TRUE
    if(missing(valType)) {
        valType <- "p.greater"
    } else {
        valType <- match.arg(valType)
    }
    typeInt <- type2int(valType)
    
    if(storage.mode(matrix)=="character")
        stop("Input must be a numeric matrix or anything that can be converted into a numeric matrix")
    
    if(storage.mode(matrix)!="double")
        storage.mode(matrix) <- "double"

    if(offset(indexList)!=0L)
        offset(indexList) <- 0L
    
    res <- .Call(C_wmw_test, matrix, indexList, typeInt)
    rownames(res) <- names(indexList)
    colnames(res) <- colnames(matrix)

    if(simplify) {
        res <- simplifyMatrix(res)
    }
    return(res)
}
#'@describeIn wmwTest \code{x} is a \code{matrix} and \code{indexList} is a \code{IndexList}
setMethod("wmwTest", c("matrix", "IndexList"),
          function(x, indexList,
                   valType, simplify = TRUE) {
              wmwTest.default(x, indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is a \code{numeric} and \code{indexList} is a \code{IndexList}
setMethod("wmwTest", c("numeric", "IndexList"),
          function(x, indexList,
                    valType, simplify = TRUE) {
              x <- matrix(x, ncol=1)
              wmwTest.default(x, indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is a \code{matrix} and \code{indexList} is a \code{GmtList}
setMethod("wmwTest", c("matrix", "GmtList"),
          function(x, indexList,
                   valType, simplify = TRUE) {
              indexList <- matchGenes(indexList, x)
              wmwTest.default(x, indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is a \code{eSet} and \code{indexList} is a \code{GmtList}
setMethod("wmwTest", c("eSet", "GmtList"),
          function(x, indexList, col="GeneSymbol",
                   valType, simplify = TRUE) {
              indexList <- matchGenes(indexList, x, col = col)
              wmwTest.default(exprs(x), indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is a \code{eSet} and \code{indexList} is a \code{numeric}
setMethod("wmwTest", c("eSet", "numeric"),
          function(x, indexList, col="GeneSymbol",
                   valType, simplify = TRUE) {
              wmwTest(exprs(x), indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is a \code{eSet} and \code{indexList} is a \code{logical}
setMethod("wmwTest", c("eSet", "logical"),
          function(x, indexList, col="GeneSymbol",
                   valType, simplify = TRUE) {
              wmwTest(exprs(x), indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is a \code{eSet} and \code{indexList} is a \code{list}
setMethod("wmwTest", c("eSet", "list"),
          function(x, indexList, col="GeneSymbol",
                   valType, simplify = TRUE) {
              indexList <- IndexList(indexList)
              wmwTest(exprs(x), indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is \code{ANY} and \code{indexList} is a \code{numeric}
setMethod("wmwTest", c("ANY", "numeric"),
          function(x, indexList, valType, simplify = TRUE) {
              indexList <- IndexList(indexList)
              wmwTest(x, indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is \code{ANY} and \code{indexList} is a \code{logical}
setMethod("wmwTest", c("ANY", "logical"),
          function(x, indexList, valType, simplify = TRUE) {
              indexList <- IndexList(indexList)
              wmwTest(x, indexList, valType=valType, simplify=simplify)
          })
#'@describeIn wmwTest \code{x} is \code{ANY} and \code{indexList} is a \code{list}
setMethod("wmwTest", c("ANY", "list"),
          function(x, indexList, valType, simplify = TRUE) {
              indexList <- IndexList(indexList)
              wmwTest(x, indexList, valType=valType, simplify=simplify)
          })


##----------------------------------------##
## wmwTestSignedGenesets
##----------------------------------------##
#'@title
#' 
#'Wilcoxon-Mann-Whitney rank sum test for high-throughput expression
#'profiling data
#'
#'@description
#'We have implemented a highly efficient Wilcoxon-Mann-Whitney rank sum
#'test for high-throughput expression profiling data. For datasets with
#'more than 100 features (genes), the function can be more than 1,000 
#'times faster than its R implementations (\code{wilcox.test} in 
#'\code{stats}, or \code{rankSumTestWithCorrelation} in \code{limma}).
#'
#'@param x A numeric matrix. All other data types (e.g. numeric vectors
#'or \code{ExpressionSet} objects) are coerced into matrix.
#'
#'@param indexList A list of integer indices (starting from 1) indicating
#'signature genes. Can be of length zero. Other data types (e.g. a list
#'of numeric or logical vectors, or a numeric or logical vector) are
#'coerced into such a list. See \code{details} below for a special case
#'using GMT files.
#'
#'@param valType The value type to be returned, allowed values
#'include \code{p.greater}, \code{p.less}, \code{abs.log10p.greater} and 
#'\code{abs.log10p.less} (one-sided tests),\code{p.two.sided}, and \code{U} 
#'statistic, and their log10 transformation variants. See details below.
#'
#'@param col a string sometimes used with a \code{eSet}
#'
#'@param simplify Logical. If not, the returning value is in matrix
#'format; if set to \code{TRUE}, the results are simplified into
#'vectors when possible (default).
#'
#'@details The basic application of the function is to test the enrichment of
#'gene sets in expression profiling data or differentially expressed
#'data (the matrix with feature/gene in rows and samples in columns).
#'
#'A special case is when \code{x} is an \code{eSet} object
#'(e.g. \code{ExpressionSet}), and \code{indexList} is a list returned
#'from \code{readGmt} function. In this case, the only requirement is
#'that one column named \code{GeneSymbol} in the \code{featureData}
#'contain gene symbols used in the GMT file. See the example below.
#'
#'Besides the conventional value types such as \sQuote{p.greater},
#'\sQuote{p.less}, \sQuote{p.two.sided} , and \sQuote{U} (the U-statistic), 
#'\code{wmwTest} (from version 0.99-1) provides further value types:
#'\code{abs.log10p.greater} and \code{log10p.less} perform log10
#'transformation on respective \emph{p}-values and give the
#'transformed value a proper sign (positive for greater than, and
#'negative for less than); \code{abs.log10p.two.sided} transforms
#'two-sided \emph{p}-values to non-negative values; and \code{Q} score
#'reports absolute log10-transformation of \emph{p}-value of the
#'two-side variant,  and gives a proper sign to it, depending on whether it is
#'rather greater than (positive) or  less than (negative). 
#'
#'@return A numeric matrix or vector containing the statistic.
#'
#'@references Barry, W.T., Nobel, A.B., and Wright, F.A. (2008). A statistical framework for testing functional namespaces in microarray data. _Annals of Applied Statistics_ 2, 286-315.
#'
#'Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test
#'accounting for inter-gene correlation. _Nucleic Acids Research_ 40(17):e133
#'
#'Zar, JH (1999). _Biostatistical Analysis 4th Edition_. Prentice-Hall International, Upper Saddle River, New Jersey.
#'
#' @import Rcpp
#' @importFrom Biobase fData exprs
#' @importClassesFrom Biobase eSet
#' 
#'@author Jitao David Zhang <jitao_david.zhang@roche.com>
#'
#'@note The function has been optimized for expression profiling data. It
#'avoids repetitive ranking of data as done by native R implementations
#'and uses efficient C code to increase the performance and control
#'memory use. Simulation studies using expression profiles of 22000
#'genes in 2000 samples and 200 gene sets suggested that the C
#'implementation can be >1000 times faster than the R
#'implementation. And it is possible to further accelerate by
#'parallel calling the function with \code{mclapply} in the \code{multicore} package.
#'
#'@seealso code{wilcox.test} in the \code{stats} package, and \code{rankSumTestWithCorrelation} in
#'the \code{limma} package.
#'
#'@examples 
#'## R-native data structures
#'set.seed(1887)
#'rd <- rnorm(1000)
#'rl <- sample(c(TRUE, FALSE), 1000, replace=TRUE)
#'wmwTest(rd, rl, valType="p.two.sided")
#'wmwTest(rd, which(rl), valType="p.two.sided")
#'rd1 <- rd + ifelse(rl, 0.5, 0)
#'wmwTest(rd1, rl, valType="p.greater")
#'wmwTest(rd1, rl, valType="U")
#'rd2 <- rd - ifelse(rl, 0.2, 0)
#'wmwTest(rd2, rl, valType="p.greater")
#'wmwTest(rd2, rl, valType="p.two.sided")
#'wmwTest(rd2, rl, valType="p.less")
#'
#'## matrix forms
#'rmat <- matrix(c(rd, rd1, rd2), ncol=3, byrow=FALSE)
#'wmwTest(rmat, rl, valType="p.two.sided")
#'wmwTest(rmat, rl, valType="p.greater")
#'
#'wmwTest(rmat, which(rl), valType="p.two.sided")
#'wmwTest(rmat, which(rl), valType="p.greater")
#'
#'## other valTypes
#'wmwTest(rmat, which(rl), valType="U")
#'wmwTest(rmat, which(rl), valType="abs.log10p.greater")
#'wmwTest(rmat, which(rl), valType="log10p.less")
#'wmwTest(rmat, which(rl), valType="abs.log10p.two.sided")
#'wmwTest(rmat, which(rl), valType="Q")
#'
#'## using ExpressionSet
#'data(sample.ExpressionSet)
#'testSet <- sample.ExpressionSet
#'fData(testSet)$GeneSymbol <- paste("GENE_",1:nrow(testSet), sep="")
#'mySig1 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.25, 0.75), replace=TRUE)
#'wmwTest(testSet, which(mySig1), valType="p.greater")
#'
#'## using integer
#'exprs(testSet)[,1L] <- exprs(testSet)[,1L] + ifelse(mySig1, 50, 0)
#'wmwTest(testSet, which(mySig1), valType="p.greater")
#'
#'## using lists
#'mySig2 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.6, 0.4), replace=TRUE)
#'wmwTest(testSet, list(first=mySig1, second=mySig2))

#'## using GMT file
#'gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
#'gmt_list <- readGmt(gmt_file)
#'
#'gss <- sample(unlist(sapply(gmt_list, function(x) x$genes)), 1000)
#'eset<-new("ExpressionSet",
#'          exprs=matrix(rnorm(10000), nrow=1000L),
#'          phenoData=new("AnnotatedDataFrame", data.frame(Sample=LETTERS[1:10])),
#'          featureData=new("AnnotatedDataFrame",data.frame(GeneSymbol=gss)))
#'esetWmwRes <- wmwTest(eset ,gmt_list, valType="p.greater")
#'summary(esetWmwRes)
#'@name wmwTest
NULL

#' @importFrom methods is
wmwTestSignedGenesets.default <- function(matrix,
                                          signedIndexList,
                                          valType=c("p.greater", "p.less", "p.two.sided", "U",
                                              "abs.log10p.greater","log10p.less","abs.log10p.two.sided",
                                              "Q"),
                                          simplify=TRUE) {
    if(!is.matrix(matrix) || !is(signedIndexList, "SignedIndexList"))
        stop("'matrix' and 'signedIndexList' must be matrix and an SignedIndexList object, respectively")
    if(missing(simplify))  simplify <- TRUE
    if(missing(valType)) {
        valType <- "p.greater"
    } else {
        valType <- match.arg(valType)
    }
    typeInt <- type2int(valType)
    
    if(storage.mode(matrix)=="character")
        stop("Input must be a numeric matrix or anything that can be converted into a numeric matrix")
    
    if(storage.mode(matrix)!="double")
        storage.mode(matrix) <- "double"

    if(offset(signedIndexList)!=0L)
        offset(signedIndexList) <- 0L
    
    res <- .Call(C_signed_wmw_test, matrix,
                 signedIndexList,
                 typeInt)
    rownames(res) <- names(signedIndexList)
    colnames(res) <- colnames(matrix)

    if(simplify) {
        res <- simplifyMatrix(res)
    }
    return(res)
}
#'@describeIn wmwTest \code{x} is a \code{matrix} and \code{indexList} is a 
#'\code{SignedIndexList}
setMethod("wmwTest", c("matrix", "SignedIndexList"), function(x, indexList, valType, simplify = TRUE) {
    wmwTestSignedGenesets.default(x, indexList, valType, simplify = simplify)
})
#'@describeIn wmwTest \code{x} is a \code{numeric} and \code{indexList} is a 
#'\code{SignedIndexList}
setMethod("wmwTest", c("numeric", "SignedIndexList"), function(x, indexList, valType, simplify = TRUE) {
              x <- matrix(x, ncol=1L)
              wmwTestSignedGenesets.default(x, indexList, valType, simplify = simplify)
})
#'@describeIn wmwTest \code{x} is a \code{eSet} and \code{indexList} is a 
#'\code{SignedIndexList}
setMethod("wmwTest", c("eSet", "SignedIndexList"), function(x, indexList, valType, simplify = TRUE) {
              wmwTestSignedGenesets.default(exprs(x), indexList, valType, simplify = simplify)
          })

##setGeneric("wmwTest",function(object, sub, alternative, statistic) standardGeneric("wmwTest"))
##setGeneric("wmwTest",function(exprs, index, alternative) standardGeneric("wmwTest"))
##setMethod("wmwTest", signature=c("ANY", "ANY", "character") , function(exprs, index, alternative) {
##  wmwTestC(exprs, index, alternative=alternative)
##})
##setMethod("wmwTest", signature=c("eSet", "gmtlist", "character"), function(exprs, index, alternative) {
##  if(!"GeneSymbol" %in% colnames(fData(exprs)))
##    stop("ExpressionSet must has 'GeneSymbol' as fData column which contains gene symbols used in the GMT files\n")
##  gb <- as.character(fData(exprs)[,"GeneSymbol"])
##  ind <- lapply(index, function(x) {
##    rind <- match(x$genes, gb)
##    rind <- rind[!is.na(rind)]
##    return(rind-1L)
##  })
##  names(ind) <- sapply(index, function(x) x$name)
##  wmwTest(exprs(exprs), ind, alternative=alternative)
##})
