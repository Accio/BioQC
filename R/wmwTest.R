## FOR REFERENCE rankSumTestWithCorrelation function from the limma package (version 3.18.13)

## authors: Gordon Smyth and Di Wu, following Zar, JD (1999) Biostatistical Analysis 4th Edition
## used under GPL(>=2) license. The function has been slihtly modified to allow reporting results

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
#' @examples
#' testNums <- 1:10
#' testSub <- rep_len(c(TRUE, FALSE), length.out=length(testNums))
#' wmwTestInR(testNums, testSub)
#' wmwTestInR(testNums, testSub, valType="p.two.sided")
#' wmwTestInR(testNums, testSub, valType="p.less")
#' wmwTestInR(testNums, testSub, valType="W")
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

valTypes <- function() names(TYPE_CODES)

type2int <- function(type) {
    if(!type %in% names(TYPE_CODES))
        stop("Should not happen! This is the wrong code")
    return(TYPE_CODES[type])
}

gmtList2IndexList.default <- function(geneSymbols, gmtList) {
    genes <- lapply(gmtList, function(x) x$genes)
    names(genes) <- sapply(gmtList, function(x) x$name)
    indList <- sapply(genes, function(g) match(g, geneSymbols))
    res <- IndexList(indList)
    return(res)
}

setMethod("gmtList2IndexList", c("character", "GmtList"), function(object, list) {
              gmtList2IndexList.default(object, list)
          })
setMethod("gmtList2IndexList", c("matrix", "GmtList"), function(object, list) {
              if(is.null(rownames(object)))
                  stop("When used to map genes in GmtList directly to rows in matrix, the matrix's row names must be gene symbols")
              symbols <- rownames(object)
              gmtList2IndexList.default(symbols, list)
          })
setMethod("gmtList2IndexList", c("eSet", "GmtList"), function(object, list, col="GeneSymbol") {
              if(!is.null(col) && !col %in% colnames(fData(object)))
                  stop("When used to map genes in GmtList directly to rows in an eSet, col must be either NULL (mapped to feature names) or a column in the fData(eset)")
              if(is.null(col)) {
                  symbols <- featureNames(object)
              } else {
                  symbols <- fData(object)[,col]
              }
              gmtList2IndexList.default(symbols, list)
          })

wmwTest.default <- function(matrix,
                            indexList,
                            valType=c("p.greater", "p.less", "p.two.sided", "U",
                                "abs.log10p.greater","log10p.less","abs.log10p.two.sided",
                                "Q"),
                            simplify=TRUE) {
    
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
    
    res <- .Call("wmw_test", matrix, indexList, typeInt)
    rownames(res) <- names(indexList)
    colnames(res) <- colnames(matrix)

    if(simplify) {
        res <- simplifyMatrix(res)
    }
    return(res)
}

setMethod("wmwTest", c("matrix", "IndexList"),
          function(object,indexList,
                   valType, simplify) {
              wmwTest.default(object, indexList, valType=valType, simplify=simplify)
          })

setMethod("wmwTest", c("numeric", "IndexList"),
          function(object, indexList,
                    valType, simplify) {
              object <- matrix(object, ncol=1)
              wmwTest.default(object, indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("matrix", "GmtList"),
          function(object, indexList,
                   valType, simplify) {
              indexList <- gmtList2IndexList(object, indexList)
              wmwTest.default(object, indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("eSet", "GmtList"),
          function(object, indexList, col="GeneSymbol",
                   valType, simplify) {
              indexList <- gmtList2IndexList(object, indexList, col=col)
              wmwTest.default(exprs(object), indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("eSet", "numeric"),
          function(object, indexList, col="GeneSymbol",
                   valType, simplify) {
              wmwTest(exprs(object), indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("eSet", "logical"),
          function(object, indexList, col="GeneSymbol",
                   valType, simplify) {
              wmwTest(exprs(object), indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("eSet", "list"),
          function(object, indexList, col="GeneSymbol",
                   valType, simplify) {
              indexList <- IndexList(indexList)
              wmwTest(exprs(object), indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("ANY", "numeric"),
          function(object, indexList, valType, simplify) {
              indexList <- IndexList(indexList)
              wmwTest(object, indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("ANY", "logical"),
          function(object, indexList, valType,simplify) {
              indexList <- IndexList(indexList)
              wmwTest(object, indexList, valType=valType, simplify=simplify)
          })
setMethod("wmwTest", c("ANY", "list"),
          function(object, indexList, valType,simplify) {
              indexList <- IndexList(indexList)
              wmwTest(object, indexList, valType=valType, simplify=simplify)
          })
##
##wmwTest <- function(x, ind.list,
##                    alternative=c("greater", "less", "two.sided", "U",
##                      "abs.log10.greater","log10.less","abs.log10.two.sided","Q"), simplify=TRUE) {
##
##    matrixObj <- formatMatrix(x)
##    matrix <- matrixObj$matrix
##    
##    indObj <- formatInd(ind.list, x, nrow(matrix))
##    indC <- indObj$indC
##    
##    typeInt <- type2int(match.arg(alternative))
##
##    res <- .Call("wmw_test", matrix, indC, typeInt)
##    
##    rownames(res) <- names(ind.list)
##    colnames(res) <- colnames(matrix)
##    
##    if(simplify) {
##        res <- simplifyMatrix(res)
##    }
##  
##  return(res)
##}

##----------------------------------------##
## wmwTestSignedGenesets
##----------------------------------------##

formatMatrix <- function(x) {
    isMatVec <- FALSE

    if(is(x, "eSet")) {
        matrix <- exprs(x)
    } else if (!is.matrix(x) & is.numeric(x)) {
        matrix <- matrix(x, ncol=1L)
        isMatVec <- TRUE
    } else if (is.matrix(x)) {
        matrix <- x
    } else {
        stop("'matrix' must be a numeric matrix, or a numeric vector, or an eSet object")
    }

    if(storage.mode(matrix)!="double")
        storage.mode(matrix) <- "double"

    return(list(matrix=matrix,
                isMatVec=isMatVec))
}

gmtlist2signedInd <- function(gmtlist, x,  posPattern="_UP$", negPattern="_DN$",
                              nomatch=c("ignore", "pos", "neg")) {
    if(is(x, "eSet")) {
        if(!"GeneSymbol" %in% colnames(fData(x)))
            stop("When ind.list is a 'gmtlist' and 'x' must is an eSet object, x's fData must contain the column 'GeneSymbol'")
        inputGenes <- fData(x)$GeneSymbol
    } else if (is.matrix(x)) {
        inputGenes <- rownames(x)
    } else {
        stop("When ind.list is a 'gmtlist', x must be either a matrix or an eSet object")
    }
    signedGenesets <- gmtlist2signedGenesets(gmtlist,
                                             posPattern=posPattern, negPattern=negPattern,
                                             nomatch=nomatch)
    indList <- sapply(signedGenesets, function(g) list(pos=match(g$pos, inputGenes),
                                                       neg=match(g$neg, inputGenes)))
    names(indList) <- names(signedGenesets)
    return(indList)
}


getCind <- function(inds, nrow) {
    if(is.logical(inds))
        inds <- which(inds)
    if(!is.numeric(inds)) {
        stop(paste("index must be either logical or numeric!\nGot ",
                   paste(head(inds), ",...", sep=" ")))
    }
    inds <- as.integer(inds[!is.na(inds)])
    if(any(inds<1L))
        stop("Indices in ind.list must be equal to or greater than one!")
    if(any(inds>nrow))
        stop("Indices out of range: they cannot exceed the row of the matrix!")
    return(inds-1L)
}

formatSignedInd <- function(signedIndList, x, nrow) {

    if(!is.list(signedIndList))
        stop("signedIndList must be a list")
    if(!all(sapply(signedIndList, function(x) length(x)==2)))
        stop("signedIndList must be a list of two items, first containing indices of positive gene sets, and the second containing negative gene sets")

    isIndVec <- length(signedIndList)==1L
    indC <- lapply(signedIndList, function(x) list(pos=getCind(x$pos, nrow=nrow),
                                                   neg=getCind(x$neg, nrow=nrow)))
    res <- list(indC=indC,
                isIndVec=isIndVec)
    return(res)
}

wmwTestSignedGenesets <- function(x, signedIndList,
                                  alternative=c("greater", "less", "two.sided", "U",
                                      "abs.log10.greater","log10.less","abs.log10.two.sided","Q"), simplify=TRUE) {
    matrixObj <- formatMatrix(x)
    matrix <- matrixObj$matrix

    indObj <- formatSignedInd(signedIndList, x, nrow(matrix))
    typeInt <- type2int(match.arg(alternative))

    res <- .Call("signed_wmw_test", matrix, posObj$indC, negObj$indC, typeInd)

    rownames(res) <- names(ind.list)
    colnames(res) <- colnames(matrix)

    if(simplify) {
        res <- simplifyWmw(res, matrixObj$isMatVec, posObj$isIndVec)
    }
    
    return(res)
}

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
