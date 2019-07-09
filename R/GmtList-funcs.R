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