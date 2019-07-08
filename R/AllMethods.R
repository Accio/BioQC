##----------------------------------------##
## standard generics
##----------------------------------------##
#' @rdname offset
#' @usage offset(object)
#' @exportMethod offset
setGeneric("offset", function(object) standardGeneric("offset"))

#'@rdname offset-set
#'@usage `offset<-`(object, value)
#'@exportMethod offset<-
setGeneric("offset<-", function(object, value) standardGeneric("offset<-"))

#' @rdname IndexList
#' @exportMethod IndexList
setGeneric("IndexList", function(object, ..., keepNA = FALSE, keepDup = FALSE, offset =1L) standardGeneric("IndexList"))

#'@rdname SignedIndexList
#'@exportMethod SignedIndexList
setGeneric("SignedIndexList", function(object, ...) standardGeneric("SignedIndexList"))

#'@rdname matchGenes
#'@exportMethod matchGenes
setGeneric("matchGenes", function(list, object, ...) standardGeneric("matchGenes"))

#'@rdname wmwTest
#'@exportMethod wmwTest
setGeneric("wmwTest", function(x, indexList, col = "GeneSymbol", valType = c("p.greater", "p.less", "p.two.sided","U","abs.log10p.greater","log10p.less","abs.log10p.two.sided","Q"), simplify = TRUE) standardGeneric("wmwTest"))


##--------------------##
## IndexList
##--------------------##
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

#' @importFrom methods new
IndexListFromList <- function(inlist, keepNA=FALSE, keepDup=FALSE, offset=1L) {
    outlist <- lapply(inlist, parseIndex, keepNA=keepNA, keepDup=keepDup)
    ## note that .Data must be the first data slot in the new command
    res <- new("IndexList", .Data=outlist, keepNA=keepNA, keepDup=keepDup, offset=as.integer(offset))
    return(res)
}

#' Convert a list to an IndexList object
#' 
#' @param object Either a list of unique integer indices, NULL and logical 
#' vectors (of same lengths), or a numerical vector or a logical vector. NA is discarded.
#' @param ...  If \code{object} isn't a list, additional vectors can go here.
#' @param keepNA Logical, whether NA indices should be kept or not. Default: FALSE (removed)
#' @param keepDup Logical, whether duplicated indices should be kept or not. Default: FALSE (removed)
#' @param offset Integer, the starting index. Default: 1 (as in the convention of R)
#' @return The function returns a list of vectors
#' @examples
#' testList <- list(GS_A=c(1,2,3,4,3),
#'                  GS_B=c(2,3,4,5),
#'                  GS_C=NULL,
#'                  GS_D=c(1,3,5,NA),
#'                  GS_E=c(2,4))
#' testIndexList <- IndexList(testList)
#' @name IndexList
#' @exportMethod IndexList
NULL

#' @rdname IndexList
setMethod("IndexList", "numeric", function(object, ..., keepNA=FALSE, keepDup=FALSE, offset=1L) {
              olist <- list(...)
              list <- c(list(object), olist)
              IndexListFromList(list, keepNA=keepNA, keepDup=keepDup, offset=offset)
          })


#' @examples
#' IndexList(c(FALSE, TRUE, TRUE), c(FALSE, FALSE, TRUE), c(TRUE, FALSE, FALSE), offset=0)
#' @rdname IndexList
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

#' @examples
#' IndexList(list(A=1:3, B=4:5, C=7:9))
#' IndexList(list(A=1:3, B=4:5, C=7:9), offset=0)
#' @rdname IndexList
setMethod("IndexList", "list", function(object, keepNA=FALSE, keepDup=FALSE, offset=1L) {
              IndexListFromList(object, keepNA=keepNA, keepDup=keepDup, offset=offset)
          })

##----------------------------------------##
## SignedIndexList
##----------------------------------------##
#' @importFrom methods new
SignedIndexListFromList <- function(inlist, keepNA=FALSE, keepDup=FALSE, offset=1L) {
  if(length(inlist)==2 && all(c("pos", "neg") %in% names(inlist))) {
    inlist <- list(inlist)
  }
  outlist <- lapply(inlist, function(x) list(pos=parseIndex(x$pos, keepNA=keepNA, keepDup=keepDup),
                                             neg=parseIndex(x$neg, keepNA=keepNA, keepDup=keepDup)))
  ## note that .Data must be the first data slot in the new command
  res <- new("SignedIndexList", .Data=outlist, keepNA=keepNA, keepDup=keepDup, offset=as.integer(offset))
  return(res)
}

#'Convert a list into a SignedIndexList
#'@name SignedIndexList
NULL
#'@param object A list of lists, each with two elements named `pos` or `neg`, can be logical vectors or integer indices
#'@param ... additional arguments, currently ignored
#'@param keepNA Logical, whether NA indices should be kept or not. Default: 
#'FALSE (removed)
#'@param keepDup Logical, whether duplicated indices should be kept or not. 
#'Default: FALSE (removed) 
#'@param offset offset; 1 if missing
#'@return A SignedIndexList, a list of lists, containing two vectors named `positive` and `negative`, 
#' which contain the indices of genes that are either positively or negatively associated with a certain
#' phenotype
#' 
#'@examples
#'myList <- list(a = list(pos = list(1, 2, 2, 4), neg = c(TRUE, FALSE, TRUE)), 
#'b = list(NA), c = list(pos = c(c(2, 3), c(1, 3))))
#'SignedIndexList(myList)
#'
#'## a special case of input is a single list with two elements, \code{pos} and \code{neg}
#'SignedIndexList(myList[[1]])
#'@rdname SignedIndexList
#'@export SignedIndexList
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
#' @name offset
#' @aliases offset,BaseIndexList-method
#' @export 
setMethod("offset", "BaseIndexList", function(object) return(object@offset))

modOffset <- function(x, diff) {
    if(is.null(x)) return(NULL)
    return(x-diff)
}

#' Set the offset of an \code{IndexList} or a \code{SignedIndexList} object
#' 
#' @param object An \code{IndexList} or a \code{SignedIndexList} object
#' @param value The value, that the offset of \code{object} is set too. If it 
#' isn't an integer, it's coerced into an integer.
#' @examples 
#' myIndexList <- IndexList(list(1:5, 2:7, 3:8), offset=1L)
#' offset(myIndexList)
#' offset(myIndexList) <- 3
#' offset(myIndexList)
#' @name offset-set
NULL

#'@rdname offset-set
setMethod("offset<-", c("IndexList", "numeric"), function(object, value) {
              value <- as.integer(value)
              diff <- object@offset - value
              object@offset <- value
              resList <- lapply(object@.Data, modOffset, diff=diff)
              object@.Data <- resList
              return(object)
          })

#'@rdname offset-set
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
## subsetting
##--------------------##
#' Subsetting GmtList object into another GmtList object
#' @param x A GmtList object
#' @param i Index to subset
#' @param drop In case only one element remains, should a list representing the single geneset returned? Default: FALSE
#' 
#' @importFrom methods new
#' 
#' @examples 
#' myGmtList <- GmtList(list(gs1=letters[1:3], gs2=letters[3:4], gs3=letters[4:5]))
#' myGmtList[1:2]
#' myGmtList[1] ## default behaviour: not dropping
#' myGmtList[1,drop=TRUE] ## force dropping
`[.GmtList` <- function(x, i, drop=FALSE) {
  if(is.character(i))
    i <- match(i, names(x))
  res <- new("GmtList", .Data=x@.Data[i])
  if(length(res)==1 && drop)
    res <- res@.Data[[1]]
  names(res) <- names(x)[i]
  return(res)
}

#' Subsetting GmtList object to fetch one gene-set
#' @param x A GmtList object
#' @param i The index to subset
#' @examples 
#' myGmtList <- GmtList(list(gs1=letters[1:3], gs2=letters[3:4], gs3=letters[4:5]))
#' myGmtList[[1]]
`[[.GmtList` <- function(x, i) {
  if(is.character(i))
    i <- match(i, names(x))
  res <- x@.Data[[i]]
  return(res)
}

##-------------------------------##
## convenience funcs for GmtList
##-------------------------------##
#' Gene-set names
#' 
#' @param x A \code{GmtList} object
#' @return Names as a vector of character strings of the same length as \code{x}
#' @export gsName
gsName <- function(x) sapply(x, function(xx) xx$name)


#' Gene-set descriptions
#' 
#' @param x A \code{GmtList} object
#' @return Descriptions as a vector of character strings of the same length as \code{x}
#' @export gsDesc
gsDesc <-  function(x) sapply(x, function(xx) xx$desc)


#' Gene-set member genes
#' 
#' @param x A \code{GmtList} object
#' @return A list of genes as character strings of the same length as \code{x}
#' @export gsGenes
gsGenes <- function(x) lapply(x, function(xx) xx$genes)


#' Gene-set gene counts
#' 
#' @param x A \code{GmtList} or similar object
#' @param uniqGenes Logical, whether only unique genes are counted
#' 
#' @return Gene counts (aka gene-set sizes) as a vector of integer of the same length as \code{x}
#' @export gsGeneCount
gsGeneCount <- function(x, uniqGenes=TRUE) {
  res <- sapply(x, function(x) {
    genes <- x$genes
    if(uniqGenes)
      genes <- unique(genes)
    return(length(genes))
  })
  return(res)
}

#' gsSize is the synonym of gsGeneCount
#' @rdname gsGeneCount
#' @export gsSize
gsSize <- function(x, uniqGenes=TRUE) gsGeneCount(x, uniqGenes=uniqGenes)


#' Whether namespace is set
#' 
#' @param x A \code{GmtList} object
#' @return Logical, whether all gene-sets have the field \code{namespace} set
#' @export hasNamespace
hasNamespace <- function(x) all(sapply(x, function(xx) !is.null(xx$namespace)))

#' Gene-set namespaces
#' 
#' @param x A \code{GmtList} object
#' @return Namespaces as a vector of character strings of the same length as \code{x}
#' @export gsNamespace
gsNamespace <- function(x) sapply(x, function(xx) xx$namespace)

#' Filter a GmtList by size
#' 
#' @param x A \code{GmtList} object
#' @param min Numeric, gene-sets with fewer genes than \code{min} will be removed
#' @param max Numeric, gene-sets with more genes than \code{max} will be removed
#' 
#' @return A \code{GmtList} object with sizes (count of genes) between \code{min} and \code{max} (inclusive).
#' @export filterBySize
filterBySize <- function(x, min, max) function(x, min, max) {
  sizes <- gsSize(x)
  isKept <- rep.int(TRUE, length(sizes))
  if(!missing(min))
    isKept[sizes<min] <- FALSE
  if(!missing(max))
    isKept[sizes>max] <- FALSE
  res <- x[isKept]
  return(res)
}

#' Set the namespace field in each gene-set within a GmtList
#' 
#' @param x A \code{GmtList} object encoding a list of gene-sets
#' @param namespace It can be either a function that applies to a \code{gene-set list} element of the object (for instance \code{function(x) x$desc} to extract description), or a vector of the same length of \code{x}, or in the special case \code{NULL}, which will erase the field namespace.
#' 
#' Note that using vectors as \code{namespace} leads to poor performance when the input object has many gene-sets.
#' 
#' @examples 
#' myGmtList <- GmtList(list(list(name="GeneSet1", desc="Namespace1", genes=LETTERS[1:3]),
#'   list(name="GeneSet2", desc="Namespace1", genes=rep(LETTERS[4:6],2)),
#'   list(name="GeneSet1", desc="Namespace1", genes=LETTERS[4:6]),
#'   list(name="GeneSet3", desc="Namespace2", genes=LETTERS[1:5])))
#' hasNamespace(myGmtList)
#' myGmtList2 <- setNamespace(myGmtList, namespace=function(x) x$desc)
#' gsNamespace(myGmtList2)
#' ## the function can provide flexible ways to encode the gene-set namespace
#' myGmtList3 <- setNamespace(myGmtList, namespace=function(x) gsub("Namespace", "C", x$desc))
#' gsNamespace(myGmtList3)
#' ## using vectors
#' myGmtList4 <- setNamespace(myGmtList, namespace=c("C1", "C1", "C1", "C2"))
#' gsNamespace(myGmtList4)
#' myGmtList2null <- setNamespace(myGmtList2, namespace=NULL)
#' hasNamespace(myGmtList2null)
#' @export
setNamespace <- function(x, namespace) {
  if(missing(namespace)) {
    stop("'namespace' must be given. It can be a function applied to each geneset-, NULL, or a vector.")
  }
  if(is.function(namespace)) {
    res <- GmtList(lapply(x, function(gs) {
      gs$namespace <- do.call(namespace, list(gs))
      return(gs)
    }))
  } else if (is.null(namespace)) {
    res <- GmtList(lapply(x, function(gs) {
      gs$namespace <- NULL
      return(gs)
    }))
  } else {
    stopifnot(length(namespace) == length(x) &&
                (is.character(namespace) || is.factor(namespace) || 
                   is.numeric(namespace) || is.logical(namespace)))
    res <- GmtList(lapply(seq(along=x), function(i) {
      gs <- x[[i]]
      gs$namespace <- namespace[i]
      return(gs)
    }))
  }
  return(res)
}

#' gsNamespace<- is the synonym of setGsNamespace
#' @rdname setGsNamespace
#' @param x A \code{GmtList} object
#' @param value \code{namespace} in \code{setGsNamespace}. It can be either a function that applies to a \code{gene-set list} element of the object (for instance \code{function(x) x$desc} to extract description), or a vector of the same length of \code{x}, or in the special case \code{NULL}, which will erase the field namespace.
#' @export gsNamespace<-
`gsNamespace<-` <- function(x, value) {
  return(setNamespace(x, namespace=value))
}

#' Set gene-set description as namespace
#' 
#' @param x A \code{GmtList} object
#' 
#' This function wrapps \code{setNamespace} to set gene-set description as namespace
#' @seealso \code{\link{setNamespace}}
#' @export
setDescAsNamespace <- function(x) {
  return(setNamespace(x, namespace = function(xx) xx$desc))
}
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

#' Show method for GmtList
#' 
#' @param object An object of the class \code{GmtList}
#' 
#' @export
setMethod("show", "GmtList", function(object) {
  str <- sprintf("A gene-set list in GMT format with %d genesets\n", length(object))

  indent <- 2
  if(hasNamespace(object)) {
    namespaces <- gsNamespace(object)
    catTbl <- table(namespaces)
    catHead <- sprintf("%s:\n",
                       ifelse(length(catTbl)>1,
                              "Namespaces", "Namespace"))
    catDetailsStr <- sprintf("%s[%d] %s (n=%d)\n",
            paste(rep(" ", indent), collapse=""),
            seq(along=catTbl),
            names(catTbl),
            unname(catTbl))
    catStr <- paste(catHead, 
                    paste(catDetailsStr, collapse=""),
                    sep="")
  } else {
    catStr <- ""
  }
  showHead <- "Gene-sets:\n"
  if(length(object)>6) {
    heads <- object[1:3]
    tails <- object[(length(object)-2):length(object)]
    shows <- c(sapply(heads, showGeneSet, indent=indent),
               paste(paste(rep(" ", indent), collapse=""), "...", sep=""),
               sapply(tails, showGeneSet, indent=indent))
  } else {
    shows <- sapply(object, showGeneSet, indent=indent)
  }
  concStr <- paste(c(str, catStr, showHead,
                     paste(shows, collapse="\n"), "\n"), collapse="")
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
        res <- sprintf("%s%s[n=%d]\n%spositive[n=%d]:%s%s\n%snegative[n=%d]:%s%s",
                       idents,
                       name,
                       posLen+negLen,
                       doubleIdents,
                       posLen,
                       paste(pos[1:pmin(posLen, nGene)],collapse=","),
                       ifelse(posLen>nGene, ",...", ""),
                       doubleIdents,
                       negLen,
                       paste(neg[1:pmin(negLen, nGene)],collapse=","),
                       ifelse(negLen>nGene, ",...", ""))
    }
    return(res)
}

#' Show method for SignedGenesets
#' 
#' @param object An object of the class \code{SignedGenesets}
#' 
#' @export
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

#' Show method for IndexList
#' 
#' @param object An object of the class \code{IndexList}
#' 
#' @export
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

##--------------------##
## show for SignedIndexList
##--------------------##
showSignedIndices <- function(indices, name, nInd=3, indent=2) {
    if(is.null(name)) name <- "NONAME"
    indices$name <- name
    showSignedGeneset(indices, nGene=nInd, indent=indent)
}

#' Show method for SignedIndexList
#' 
#' @param object An object of the class \code{SignedIndexList}
#' 
#' @export
setMethod("show", "SignedIndexList", function(object) {
              str <- sprintf("A list of %d signed indices with offset=%d", length(object), object@offset)
              opts <- sprintf("Options: NA removed: %s; duplicates removed: %s", !object@keepNA, !object@keepDup)
              indent <- 2
              if(length(object)>6) {
                  heads <- object[1:3]
                  tails <- object[(length(object)-2):length(object)]
                  shows <- c(sapply(seq(along=heads), 
                                    function(i) showSignedIndices(heads[[i]], names(heads)[i], indent=indent)),
                             paste(paste(rep(" ", indent), collapse=""), "...", sep=""),
                             sapply(seq(along=tails), function(i) 
                               showSignedIndices(tails[[i]], names(tails)[i], indent=indent)))
              } else {
                  shows <- sapply(seq(along=object), function(i) 
                    showSignedIndices(object[[i]], names(object)[i], indent=indent))
              }
              concStr <- paste(c(str, opts, shows, ""), collapse="\n")
              cat(concStr)
          })
