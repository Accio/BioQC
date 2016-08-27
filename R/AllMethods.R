##----------------------------------------##
## standard generics
##----------------------------------------##
setGeneric("offset", function(object) standardGeneric("offset"))
setGeneric("offset<-", function(object, value) standardGeneric("offset<-"))
setGeneric("IndexList", function(object, ...) standardGeneric("IndexList"))
setGeneric("SignedIndexList", function(object, ...) standardGeneric("SignedIndexList"))
setGeneric("matchGenes", function(list, object, ...) standardGeneric("matchGenes"))
setGeneric("wmwTest", function(object, indexList, ...) standardGeneric("wmwTest"))


##--------------------##
## IndexList
##--------------------##
#' Convert a list to an IndexList object
#' 
#' @param list A list of unique integer indices or NULL, or logical vectors of same lengths. NA is discarded. 
#' @param keepNA Logical, whether NA indices should be kept or not. Default: FALSE (removed)
#' @param keepDup Logical, whether duplicated indices should be kept or not. Default: FALSE (removed)
#' @param offset Integer, the starting index. Default: 1 (as in the convention of R)
#' 
#' @examples
#' testList <- list(GS_A=c(1,2,3,4,3),
#'                  GS_B=c(2,3,4,5),
#'                  GS_C=NULL,
#'                  GS_D=c(1,3,5,NA),
#'                  GS_E=c(2,4))
#' testIndexList <- IndexList(testList)

parseIndex <- function(x, keepNA=FALSE, keepDup=FALSE) {
    if(is.null(x))
        return(x)
    if(is.logical(x))
        x <- which(x)
    x <- as.integer(x)
    if(!keepNA)
        x <- x[!is.na(x)]
    if(!keepDup)
        x <- unique(x)
    return(x)
}
IndexListFromList <- function(inlist, keepNA=FALSE, keepDup=FALSE, offset=1L) {
    outlist <- lapply(inlist, parseIndex, keepNA=keepNA, keepDup=keepDup)
    ## note that .Data must be the first data slot in the new command
    res <- new("IndexList", .Data=outlist, keepNA=keepNA, keepDup=keepDup, offset=as.integer(offset))
    return(res)
}


#' Convert several numeric vectors into an index list
#'
#' @param object A integer (numeric) vector
#' @param ... Other integer (numeric) vector(s)
#' @param offset offset; 1 if missing
#'
#' @examples
#' IndexList(1:3, 4:5, 7:9)
#' IndexList(1:3, 4:5, 7:9, offset=0)
setMethod("IndexList", "numeric", function(object, ..., keepNA=FALSE, keepDup=FALSE, offset=1L) {
              olist <- list(...)
              list <- c(list(object), olist)
              IndexListFromList(list, keepNA=keepNA, keepDup=keepDup, offset=offset)
          })

#' Convert several logical vectors into an index list
#'
#' @param object A logical vector
#' @param ... Other logical vector(s) of the same length as \code{object}
#' @param offset offset; 1 if missing
#'
#' @examples
#' IndexList(c(FALSE, TRUE, TRUE), c(FALSE, FALSE, TRUE), c(TRUE, FALSE, FALSE), offset=0)
setMethod("IndexList", "logical", function(object, ..., keepNA=FALSE, keepDup=FALSE, offset=1L) {
              olist <- list(...)
              list <- c(list(object), olist)
              if(!all(sapply(list, is.logical)))
                  stop("all input must be logical vectors")
              if(!length(unique(length(sapply(list, length))))==1)
                  stop("all input must be of the same length")
              list <- lapply(list, which)
              IndexListFromList(list, keepNA=keepNA, keepDup=keepDup, offset=offset)
          })

#' Convert a list of numeric (integer) vectors into an index list
#'
#' @param object A list of numeric (integer) vectors
#' @param offset offset; 1 if missing
#'
#' @examples
#' IndexList(list(A=1:3, B=4:5, C=7:9))
#' IndexList(list(A=1:3, B=4:5, C=7:9), offset=0)
setMethod("IndexList", "list", function(object, keepNA=FALSE, keepDup=FALSE, offset=1L) {
              IndexListFromList(object, keepNA=keepNA, keepDup=keepDup, offset=offset)
          })

##----------------------------------------##
## SignedIndexList
##----------------------------------------##
SignedIndexListFromList <- function(inlist, keepNA=FALSE, keepDup=FALSE, offset=1L) {
    outlist <- lapply(inlist, function(x) list(pos=parseIndex(x$pos, keepNA=keepNA, keepDup=keepDup),
                                               neg=parseIndex(x$neg, keepNA=keepNA, keepDup=keepDup)))
    ## note that .Data must be the first data slot in the new command
    res <- new("SignedIndexList", .Data=outlist, keepNA=keepNA, keepDup=keepDup, offset=as.integer(offset))
    return(res)
}

setMethod("SignedIndexList", "list", function(object, keepNA=FALSE, keepDup=FALSE, offset=1L) {
     SignedIndexListFromList(object, keepNA=keepNA, keepDup=keepDup, offset=offset)
})

##--------------------##
## offset
##--------------------##
#' Get offset from an IndexList object
#'
#' @param object An IndexList object
#'
#' @examples
#' myIndexList <- IndexList(list(1:5, 2:7, 3:8), offset=1L)
#' offset(myIndexList)
setMethod("offset", "BaseIndexList", function(object) return(object@offset))
modOffset <- function(x, diff) {
    if(is.null(x)) return(NULL)
    return(x-diff)
}
setMethod("offset<-", c("IndexList", "numeric"), function(object, value) {
              value <- as.integer(value)
              diff <- object@offset - value
              object@offset <- value
              resList <- lapply(object@.Data, modOffset, diff=diff)
              object@.Data <- resList
              return(object)
          })
setMethod("offset<-", c("SignedIndexList", "numeric"), function(object, value) {
    value <- as.integer(value)
    diff <- object@offset - value
    object@offset <- value
    resList <- lapply(object@.Data, function(x) list(pos=modOffset(x$pos, diff),
                                                     neg=modOffset(x$neg, diff)))
    object@.Data <- resList
    return(object)
})
##--------------------##
## show for GmtList
##--------------------##
showGeneSet <- function(geneset, nGene=3, indent=2) {
    genes <- geneset$genes
    geneLen <- length(genes)
    sprintf("%s%s (%sn=%d): %s%s",
            paste(rep(" ", indent), collapse=""),
            geneset$name,
            ifelse(is.null(geneset$desc), "", paste(geneset$desc, ",", sep="")),
            geneLen,
            paste(geneset$genes[1:pmin(nGene, geneLen)],collapse=","),
            ifelse(geneLen>nGene, ",...", ""))
}
setMethod("show", "GmtList", function(object) {
              str <- sprintf("A gene-set list in gmt format with %d genesets", length(object))
              indent <- 2
              if(length(object)>6) {
                  heads <- object[1:3]
                  tails <- object[(length(object)-2):length(object)]
                  shows <- c(sapply(heads, showGeneSet, indent=indent),
                             paste(paste(rep(" ", indent), collapse=""), "...", sep=""),
                             sapply(tails, showGeneSet, indent=indent))
              } else {
                  shows <- sapply(object, showGeneSet, indent=indent)
              }
              concStr <- paste(c(str, shows, ""), collapse="\n")
              cat(concStr)
          })


##--------------------##
## show for SignedGenesets
##--------------------##
showSignedGeneset <- function(geneset, nGene=3, indent=2) {
    name <- geneset$name
    pos <- geneset$pos
    neg <- geneset$neg
    posLen <- length(pos)
    negLen <- length(neg)
    idents <- paste(rep(" ", indent), collapse="")
    doubleIdents <- paste(rep(" ", indent*2), collapse="")
    if(posLen+negLen==0) {
        res <- sprintf("%s%s (no genes)", idents, name)
    } else {
        res <- sprintf("%s%s[n=%d]\n%spositive[n=%d]:%s\n%snegative[n=%d]:%s",
                        idents,
                        name,
                        posLen+negLen,
                        doubleIdents,
                        posLen, 
                        paste(pos[1:pmin(posLen, nGene)],collapse=","),
                        doubleIdents,
                        negLen,
                        paste(neg[1:pmin(negLen, nGene)],collapse=","))
    }
    return(res)
}

setMethod("show", "SignedGenesets", function(object) {
              str <- sprintf("A list of %d signed gene-sets", length(object))
              indent <- 2
              if(length(object)>6) {
                  heads <- object[1:3]
                  tails <- object[(length(object)-2):length(object)]
                  shows <- c(sapply(heads, showSignedGeneset, indent=indent),
                             paste(paste(rep(" ", indent), collapse=""), "...", sep=""),
                             sapply(tails, showSignedGeneset, indent=indent))
              } else {
                  shows <- sapply(object, showSignedGeneset, indent=indent)
              }
              concStr <- paste(c(str, shows, ""), collapse="\n")
              cat(concStr)
          })

##--------------------##
## show for IndexList
##--------------------##
showIndices <- function(indices, name, nInd=3, indent=2) {
    if(is.null(name))
        name <- "NONAME"
    len <- length(indices)
    if(is.null(indices)) {
        isNA <- FALSE
    } else {
        isNA <- is.na(indices)
    }
    sprintf("%s%s (n=%d%s): %s%s",
            paste(rep(" ", indent), collapse=""),
            name,
            len,
            ifelse(any(isNA), sprintf(", with %d NA", sum(isNA)), ""),
            paste(indices[1:pmin(len, nInd)],collapse=","),
            ifelse(len>nInd, ",...", ""))
}
setMethod("show", "IndexList", function(object) {
              str <- sprintf("A list of %d indices with offset=%d", length(object), object@offset)
              opts <- sprintf("Options: NA removed: %s; duplicates removed: %s", !object@keepNA, !object@keepDup)
              indent <- 2
              if(length(object)>6) {
                  heads <- object[1:3]
                  tails <- object[(length(object)-2):length(object)]
                  shows <- c(sapply(seq(along=heads), function(i) showIndices(heads[[i]], names(heads)[i], indent=indent)),
                             paste(paste(rep(" ", indent), collapse=""), "...", sep=""),
                             sapply(seq(along=tails), function(i) showIndices(tails[[i]], names(tails)[i], indent=indent)))
              } else {
                  shows <- sapply(seq(along=object), function(i) showIndices(object[[i]], names(object)[i], indent=indent))
              }
              concStr <- paste(c(str, opts, shows, ""), collapse="\n")
              cat(concStr)
          })
