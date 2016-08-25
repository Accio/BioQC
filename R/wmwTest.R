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
#' @param alternative two.sided, less, or greater
#' @param statistic Logical, if \code{TRUE}, the U-statistic of the Wilcoxon-Mann-Whitney test is returned; otherwise the p-value is returned
#' 
#' @examples
#' testNums <- 1:10
#' testSub <- rep_len(c(TRUE, FALSE), length.out=length(testNums))
#' wmwTestInR(testNums, testSub)
wmwTestInR <- function(x, sub, alternative=c("two.sided", "less", "greater"), statistic=FALSE) {
    if(is.numeric(sub)) {
        tmp <- rep(FALSE, length(x))
        tmp[sub] <- TRUE
        sub <- tmp
    }
    if(!is.logical(sub))
        stop("sub must be either numeric indices or logical")
    if(!any(sub)) return(ifelse(statistic, 0, 1))
    wt <- wilcox.test(x[sub],
                      x[!sub], alternative=alternative, exact=FALSE)
    return(ifelse(statistic, wt$statistic, wt$p.value))
}

## type2int and formatMatrixInd are helper functions for wmwTest and wmwSignedTest
TYPE_CODES <- c("greater"=0L, "less"=1L,
                "two.sided"=2L, "U"=3L,
                "abs.log10.greater"=4L,
                "log10.less"=5L,
                "abs.log10.two.sided"=6L,
                "Q"=7L)

type2int <- function(type) {
    if(!type %in% names(TYPE_CODES))
        stop("Should not happen! This is the wrong code")
    return(TYPE_CODES[type])
}

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

gmtlist2ind <- function(ind.list, x) {
    if(is(x, "eSet")) {
        if(!"GeneSymbol" %in% colnames(fData(x)))
            stop("When ind.list is a 'GmtList' and 'x' must is an eSet object, x's fData must contain the column 'GeneSymbol'")
        inputGenes <- fData(x)$GeneSymbol
    } else if (is.matrix(x)) {
        inputGenes <- rownames(x)
    } else {
        stop("When ind.list is a 'GmtList', x must be either a matrix or an eSet object")
    }
    genes <- lapply(ind.list, function(x) x$genes)
    names(genes) <- sapply(ind.list, function(x) x$name)
    indList <- sapply(genes, function(g) match(g, inputGenes))
    return(indList)
}
    

formatInd <- function(ind.list, x, nrow) {
    isIndVec <- FALSE

    if(is.numeric(ind.list) || is.logical(ind.list)) {
        ind.list <- list(ind.list)
        isIndVec <- TRUE
    } else if (is(ind.list, "GmtList")) {
        ind.list <- gmtlist2ind(ind.list, x)
        isInVec <- length(ind.list)==1
    }

    indC <- lapply(ind.list, function(x) getCind(x, nrow=nrow))

    res <- list(indC=indC,
                isIndVec=isIndVec)
    return(res)
}

simplifyMatrix <- function(matrix) {
    if(nrow(matrix)==1L & ncol(matrix)==1L)  {
        matrix <- matrix[1L,1L]
    } else if (nrow(matrix)==1L) {
        matrix <- matrix[1L,]
    } else if (ncol(matrix)==1L)  {
        matrix <- matrix[,1L]
    } 
    return(matrix)
}

wmwTest.default <- function(matrix,
                            indexListInC,
                            valType=c("greater", "less", "two.sided", "U",
                                "abs.log10.greater","log10.less","abs.log10.two.sided",
                                "Q"),
                            simplify=TRUE) {

    if(storage.matrix(matrix)=="character")
        stop("Input must be a numeric matrix or anything that can be converted into a numeric matrix")
    
    if(!storage.mode(matrix)!="double")
        storage.matrix(matrix) <- "double"
    
    typeInt <- type2int(match.arg(valType))
    res <- .Call("wmw_test", matrix, indexListInC, typeInt)
    rownames(res) <- names(indexListInC)
    colnames(res) <- colnames(matrix)

    if(simplify) {
        res <- simplifyMatrix(res)
    }
    return(res)
}
                            
                            
wmwTest <- function(x, ind.list,
                    alternative=c("greater", "less", "two.sided", "U",
                      "abs.log10.greater","log10.less","abs.log10.two.sided","Q"), simplify=TRUE) {

    matrixObj <- formatMatrix(x)
    matrix <- matrixObj$matrix
    
    indObj <- formatInd(ind.list, x, nrow(matrix))
    indC <- indObj$indC
    
    typeInt <- type2int(match.arg(alternative))

    res <- .Call("wmw_test", matrix, indC, typeInt)
    
    rownames(res) <- names(ind.list)
    colnames(res) <- colnames(matrix)
    
    if(simplify) {
        res <- simplifyMatrix(res)
    }
  
  return(res)
}

##----------------------------------------##
## wmwTestSignedGenesets
##----------------------------------------##

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

formatSignedInd <- function(signedIndList, x, nrow) {

    if(!is.list(signedIndList))
        stop("signedIndList must be a list")
    if(!all(sapply(signedIndList, function(x) length(x)==2)))
        stop("signedIndList must be a list of two items, first containing indices of positive gene sets, and the second containing negative gene sets")

    isIndVec <- legnth(signedIndList)==1L
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
