## validity functions
isValidGmtList <- function(object) {
    isList <- sapply(object, is.list)
    if(!all(isList))
        return(paste("Not all items in the list is a list. \n",
                     "The following items are not:\n",
                     paste(which(!isList), collapse=","),"\n", sep=""))
    isGmt <- sapply(object, function(x) {
                        hasNames <- all(c("name", "desc", "genes") %in% names(x))
                        isValidName <- !is.null(x$name) & length(x$name)==1
                        isValidDesc <- is.null(x$desc) || length(x$desc)==1
                        return(hasNames & isValidName & isValidDesc)
                    })
    if(!all(isGmt))
        return(paste("Not all items in the list have three mandatory fields 'name', 'desc', and 'genes'",
                     "Following items are not:\n",
                     paste(which(!isGmt), collapse=","), "\n", sep=""))
    return(TRUE)
}


isValidSignedGenesets<- function(object) {
    isList <- sapply(object, is.list)
    if(!all(isList))
        return(paste("Not all items in the list is a list. \n",
                     "The following items are not:\n",
                     paste(which(!isList), collapse=","),"\n", sep=""))
    isSigned <- sapply(object, function(x) all(c("name", "pos", "neg") %in% names(x)))
    if(!all(isSigned))
        return(paste("Not all items in the list have three mandatory fields: 'name', 'pos' and 'neg'\n",
                     "Following items are not:\n",
                     paste(which(!isSigned), collapse=","), "\n", sep=""))
    isPosChar <- sapply(object, function(x) is.null(x$pos) || is.character(x$pos))
    isNegChar <- sapply(object, function(x) is.null(x$neg) || is.character(x$neg))
    if(!all(isPosChar & isNegChar))
        return(paste("Not all items in the list have character vectors in the fields 'pos' and 'neg'\n",
                     "Following items are not:\n",
                     paste(which(!(isPosChar & isNegChar)), collapse=","), "\n", sep=""))
    return(TRUE)
}

isValidIndexList <- function(object) {
    if(!(length(object@offset)==1L && is.integer(object@offset)))
        return(sprintf("offset must be a single integer. Its value now:%d; class:%s; length:%d",
                       object@offset,
                       class(object@offset),
                       length(object@offset)))

    isList <- is.list(object)
    if(!isList)
        return(paste("object must be a list of indices starting from 1"))
    isInd <- sapply(object, function(x) is.null(x) || is.integer(x) & all(x)>=1)
    if(!all(isInd))
        return(paste("object must be a list of indices starting from 1\n",
                     sprintf("Followings are not: %s", which(!isInd))))
    return(TRUE)
}

## Class definitions
setClass("GmtList", contains="list", validity=isValidGmtList)
setClass("SignedGenesets", contains="list", validity=isValidSignedGenesets)
setClass("IndexList",
         representation=list("offset"="integer"),
         prototype=prototype(offset=1L),
         contains="list", validity=isValidIndexList)

## Constructors

#' Convert a list to a GmtList object
#' @param list A list of genesets; each geneset is a list of at least three fields: 'name', 'desc', and 'genes'. 'name' and 'desc' contains one character string ('desc' can be NULL while 'name' cannot), and 'genes' can be either NULL or a character vector.
#'
#' @seealso If a list of gene symbols need to be converted into a GmtList, use 'as.gmtlist' instead
#' 
#' @examples
#' testList <- list(list(name="GS_A", desc=NULL, genes=LETTERS[1:3]),
#'                  list(name="GS_B", desc="gene set B", genes=LETTERS[1:5]),
#'                  list(name="GS_C", desc="gene set C", genes=NULL))
#' testGmt <- GmtList(testList)
GmtList <- function(list) {
    res <- new("GmtList", .Data=list)
    return(res)
}

#' Convert a list to a SignedGenesets object
#' @param list A list of genesets; each geneset is a list of at least three fields: 'name', 'pos', and 'neg'. 'name' contains one non-null character string, and both 'pos' and 'neg' can be either NULL or a character vector.
#'
#' @seealso \code{GmtList}
#' 
#' @examples
#' testList <- list(list(name="GS_A", pos=NULL, neg=LETTERS[1:3]),
#'                  list(name="GS_B", pos=LETTERS[1:5], neg=LETTERS[7:9]),
#'                  list(name="GS_C", pos=LETTERS[1:5], neg=NULL),
#'                  list(name="GS_D", pos=NULL, neg=NULL))
#' testSigndGS <- SignedGenesets(testList)
SignedGenesets <- function(list) {
    res <- new("SignedGenesets", .Data=list)
    return(res)
}

#' Convert a list to an IndexList object
#' 
#' @param list A list of indices in integer (or numeric) or NULL. 
#' 
#' @examples
#' testList <- list(GS_A=c(1,2,3,4),
#'                  GS_B=c(2,3,4,5),
#'                  GS_C=NULL,
#'                  GS_D=c(1,3,5,NA),
#'                  GS_E=c(2,4))
#' testIndexList <- IndexList(testList)
IndexListFromList <- function(list, offset=1L) {
    list <- lapply(list, function(x) {
                       if(is.null(x))
                           return(x)
                       x <- as.integer(x)
                       return(x)
                   })
    res <- new("IndexList", offset=as.integer(offset))
    res@.Data <- list
    return(res)
}

