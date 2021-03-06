#' Function to validate a GmtList object
#' @param object A GmtList object
#' Use \code{setValidity("GmtList", "isValidGmtList")} to check integrity of GmtList objects.
#' It can be very slow, therefore the feature is not turned on by default
#' @export
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

#' Function to validate a SignedGenesets object
#' @param object A SignedGenesets object
#' Use \code{setValidity("SignedGenesets", "isValidSignedGenesets")} to check integrity of SignedGenesets objects.
#' It can be very slow, therefore the feature is not turned on by default
#' @export
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

#' Function to validate a BaseIndexList  object
#' @param object A BaseIndexList  object
#' Use \code{setValidity("BaseIndexList", "isValidBaseIndexList")} to check integrity of BaseIndexList  objects.
#' It can be very slow, therefore the feature is not turned on by default
#' @export
isValidBaseIndexList <- function(object) {
  if(!(length(object@offset)==1L && is.integer(object@offset)))
    return(sprintf("offset must be a single integer. Its value now:%d; class:%s; length:%d",
                   object@offset,
                   class(object@offset),
                   length(object@offset)))
  
  isList <- is.list(object)
  if(!isList)
    return(paste("object must be a list of indices starting from 1"))
  return(TRUE)
}

#' Function to validate an IndexList object
#' @param object an IndexList object
#' Use \code{setValidity("BaseIndexList", "isValidBaseIndexList")} to check integrity of IndexList objects.
#' It can be very slow, therefore the feature is not turned on by default
#' @export
isValidIndexList <- function(object) {
  isInd <- sapply(object, function(x) is.null(x) || is.integer(x))
  if(!all(isInd))
    return(paste("object must be a list of indices starting from 1\n",
                 sprintf("Followings are not: %s", which(!isInd))))
  return(TRUE)
}

#' Function to validate a SignedIndexList object
#' @param object a SignedIndexList object
#' Use \code{setValidity("SignedIndexList", "isValidSignedIndexList")} to check integrity of SignedIndexList objects.
#' It can be very slow, therefore the feature is not turned on by default
#' @export
isValidSignedIndexList <- function(object) {
  isSigned <- sapply(object, function(x)
    all(c("pos", "neg") %in% names(x)))
  if(!all(isSigned))
    return(paste("In some items there are no pos/neg fields\n",
                 sprintf("Following are not: %s", which(!isSigned))))
  isPosInd <- sapply(object, function(x) is.null(x$pos) || is.integer(x$pos))
  isNegInd <- sapply(object, function(x) is.null(x$neg) || is.integer(x$neg))
  if(!all(isPosInd) & !all(isNegInd))
    return("Indices are not valid in following cses\n")
  return(TRUE)
}

## Class definitions

#' An S4 class to hold geneset in the GMT file in a list, each item in the list 
#' is in in turn a list containing following items: name, desc, and genes.
#' @exportClass GmtList
setClass("GmtList", 
         contains="list")

#' An S4 class to hold signed genesets, each item in the list is in in turn a 
#' list containing following items: name, pos, and neg.
#' @exportClass SignedGenesets
setClass("SignedGenesets", contains="list")


#' An S4 class to hold a list of indices, with the possibility to specify the 
#' offset of the indices. IndexList and SignedIndexList extend this class
#'
#' @slot offset An integer specifying the value of first element. Default 1
#' @slot keepNA Logical, whether NA is kept during construction
#' @slot keepDup Logical, whether duplicated values are kept during construction

setClass("BaseIndexList",
         representation=list("offset"="integer",
             "keepNA"="logical",
             "keepDup"="logical"),
         prototype=prototype(offset=1L, keepNA=FALSE, keepDup=FALSE),
         contains="list")

#' An S4 class to hold a list of integers as indices, with the possibility to specify the offset of the indices
#'
#' @slot offset An integer specifying the value of first element. Default 1
#' @slot keepNA Logical, whether NA is kept during construction
#' @slot keepDup Logical, whether duplicated values are kept during construction
#' @name IndexList-class
#' @exportClass IndexList
setClass("IndexList", contains="BaseIndexList")

#'An S4 class to hold a list of signed integers as indices, with the possibility to specify the offset of the indices
#' @slot offset An integer specifying the value of first element. Default 1
#' @slot keepNA Logical, whether NA is kept during construction
#' @slot keepDup Logical, whether duplicated values are kept during construction
#' @exportClass SignedIndexList
setClass("SignedIndexList", contains="BaseIndexList")

##----------------------------------------##
## set validity checking functions
##----------------------------------------##

## they are disabled by default, since running them over large lists is slow
## setValidity("GmtList", isValidGmtList)
## setValidity("SignedGenesets", isValidSignedGenesets)
## setValidity("BaseIndexList", isValidBaseIndexList)
## setValidity("IndexList", isValidIndexList)
## setValidity("SignedIndexList", isValidSignedIndexList)

##----------------------------------------##
## Constructors
##----------------------------------------##

#' Convert a list to a GmtList object
#' @param list A list of genesets; each geneset is a list of at least three fields: 'name', 'desc', and 'genes'. 'name' and 'desc' contains one character string ('desc' can be NULL while 'name' cannot), and 'genes' can be either NULL or a character vector. In addition, 'namespace' is accepted to represent the namespace.
#'
#' For convenience, the function also accepts a list of character vectors, each containing a geneset. In this case, the function works as a wrapper of \code{as.GmtList}
#' 
#' @seealso If a list of gene symbols need to be converted into a GmtList, use 'as.GmtList' instead
#' 
#' @examples
#' testList <- list(list(name="GS_A", desc=NULL, genes=LETTERS[1:3]),
#'                  list(name="GS_B", desc="gene set B", genes=LETTERS[1:5]),
#'                  list(name="GS_C", desc="gene set C", genes=NULL))
#' testGmt <- GmtList(testList)
#'
#' # as wrapper of as.GmtList
#' testGeneList <- list(GS_A=LETTERS[1:3], GS_B=LETTERS[1:5], GS_C=NULL)
#' testGeneGmt <- GmtList(testGeneList)
#' 
#' @importFrom methods new
#' @export
GmtList <- function(list) {
    isGeneSymbols <- all(sapply(list, function(x) is.null(x) || is.character(x)))
    if(isGeneSymbols) {
        return(as.GmtList(list))
    } else {
        res <- new("GmtList", .Data=list)
        names(res) <- names(list)
        return(res)
    }
}

#' Convert a list to a SignedGenesets object
#' @param list A list of genesets; each geneset is a list of at least three fields: 'name', 'pos', and 'neg'. 'name' contains one non-null character string, and both 'pos' and 'neg' can be either NULL or a character vector.
#'
#' @seealso \code{GmtList}
#' 
#' @importFrom methods new
#' @examples
#' testList <- list(list(name="GS_A", pos=NULL, neg=LETTERS[1:3]),
#'                  list(name="GS_B", pos=LETTERS[1:5], neg=LETTERS[7:9]),
#'                  list(name="GS_C", pos=LETTERS[1:5], neg=NULL),
#'                  list(name="GS_D", pos=NULL, neg=NULL))
#' testSigndGS <- SignedGenesets(testList)
#' 
#' @export
SignedGenesets <- function(list) {
    res <- new("SignedGenesets", .Data=list)
    return(res)
}

