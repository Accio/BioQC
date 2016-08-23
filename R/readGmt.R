##setClass("gsitem",
##         representation(name="character",
##                        desc="character",
##                        genes="character",
##                        weights="numeric",
##                        indices="integer"))
##setClass("gslist",  contains="list")
##setAs(from="list", to="gslist", def=function(from, to) {
##    new(to, from)
##})
##as.data.frame.gslist <- function(x, row.names, optional,...) {
##    genes <- gsGenes(x)
##    sets <- rep(gsName(x), sapply(genes, length))
##    unique(data.frame(GeneSet=unlist(sets), Genes=unlist(genes)))
##}
##setAs(from="gslist", to="data.frame", def=function(from,to) {
##    return(as.data.frame.gslist(from))
##})
##setAs(from="data.frame", to="gslist", def=function(from, to) {
##    ## TODO: implement from data.frame to gslist
##    stop("Not implemented yet.")
##})
##as.matrix.gslist <- function(x,...) { ## r/c: genes/gene sets
##    genes <- gsGenes(x)
##    uniqGenes <- gsUniqGenes(x)
##    res <- sapply(genes, function(xx) uniqGenes %in% xx)
##    dimnames(res) <- list(Genes=uniqGenes, GeneSets=gsName(x))
##    return(res)
##}
##setAs(from="gslist", to="matrix", def=function(from,to) {
##    return(as.matrix(from))
##})
##gsApply <- function(X, FUN, ...) {
##    res <- lapply(X, FUN, ...)
##    return(as(res, "gslist"))
##}
##tfidf <- function(x) {
##    mat <- as.matrix(x)
##    denom <- apply(mat, 2L, max)
##    tf <- 0.5+0.5*t(t(mat)/denom)
##
##    tappear <- apply(mat, 1L, function(x) sum(x!=0, na.rm=TRUE))
##    idf <- log(ncol(mat)/tappear)
##
##    return(tf*idf)
##}
##gsTfIdf <- function(x) {
##    tfidfMat <- tfidf(x)
##    genes <- gsGenes(x)
##    weights <- lapply(seq(along=genes), function(i) {
##        tfidfMat[genes[[i]],i]
##    })
##    gsWeights(x) <- weights
##    return(x)
##    
##}
##setGeneric("gsName", function(obj) standardGeneric("gsName"))
##setGeneric("gsDesc", function(obj) standardGeneric("gsDesc"))
##setGeneric("gsGenes", function(obj) standardGeneric("gsGenes"))
##setGeneric("gsWeights", function(obj) standardGeneric("gsWeights"))
##setGeneric("gsClearWeights", function(obj) standardGeneric("gsClearWeights"))
##setGeneric("gsHasWeights", function(obj) standardGeneric("gsHasWeights"))
##setGeneric("gsIndices", function(obj) standardGeneric("gsIndices"))
##setMethod("gsName", "gsitem", function(obj) obj@name)
##setMethod("gsDesc", "gsitem", function(obj) obj@desc)
##setMethod("gsGenes", "gsitem", function(obj) obj@genes)
##setMethod("gsWeights", "gsitem", function(obj) obj@weights)
##setMethod("gsHasWeights", "gsitem", function(obj) !all(obj@weights==1L))
##
##setMethod("gsIndices", "gsitem", function(obj) obj@indices)
##setMethod("length", "gsitem", function(x) length(x@genes)) 
##
##setMethod("gsName", "gslist", function(obj) sapply(obj, gsName))
##setMethod("gsDesc", "gslist", function(obj) {
##    res <- sapply(obj, gsDesc)
##    names(res) <- gsName(obj)
##    return(res)
##})
##setMethod("gsGenes", "gslist", function(obj) {
##    res <- lapply(obj, gsGenes)
##    names(res) <- gsName(obj)
##    return(res)
##})
##setMethod("gsWeights", "gslist", function(obj) {
##    res <- lapply(obj, gsWeights)
##    names(res) <- gsName(obj)
##    return(res)
##})
##setMethod("gsHasWeights", "gslist", function(obj) {
##    res <- sapply(obj, gsHasWeights)
##    names(res) <- gsName(obj)
##    return(res)
##})
##setMethod("gsIndices", "gslist", function(obj) {
##    res <- sapply(obj, gsIndices)
##    names(res) <- gsName(obj)
##    return(res)
##})
##          
##          
##setGeneric("gsName<-", function(obj,value) standardGeneric("gsName<-"))
##setGeneric("gsDesc<-", function(obj,value) standardGeneric("gsDesc<-"))
##setGeneric("gsGenes<-", function(obj,value) standardGeneric("gsGenes<-"))
##setGeneric("gsWeights<-", function(obj,value) standardGeneric("gsWeights<-"))
##setGeneric("gsIndices<-", function(obj,value) standardGeneric("gsIndices<-"))
##setMethod("gsName<-", c("gsitem","character"), function(obj,value) {
##    obj@name <- value[1]
##    return(obj)
##})
##setMethod("gsDesc<-", c("gsitem", "character"), function(obj,value) {
##    obj@desc <- value[1]
##    return(obj)
##})
##setMethod("gsGenes<-", c("gsitem","character"), function(obj,value){
##    obj@genes <- unique(value)
##    return(obj)
##})
##setMethod("gsWeights<-", c("gsitem","numeric"), function(obj,value) {
##    if(length(value)!=length(obj) && length(value)!=1)
##        stop("weights must be of the same length as the genes\n")
##    value <- rep(value, length.out=length(obj))
##    obj@weights <- unname(value)
##    return(obj)
##})
##setMethod("gsWeights<-", c("gslist", "list"), function(obj, value) {
##    if(length(value)!=length(obj))
##        stop("weights must be a list of (numeric) weights of the same length as the gslist\n")
##    gsApply(seq(along=value), function(i) {
##        gsWeights(obj[[i]]) <- value[[i]]
##        return(obj[[i]])
##    })
##})
##setMethod("gsIndices<-", c("gsitem","integer"), function(obj,value) {
##    if(length(value)!=length(obj))
##        stop("indices must be of the same length as the genes\n")
##    obj@indices <- value
##    return(obj)
##})
##setMethod("gsClearWeights", "gsitem", function(obj) {
##    gsWeights(obj) <- 1L
##    return(obj)
##})
##setMethod("gsClearWeights", "gslist", function(obj) {
##    gsApply(obj, function(x) gsClearWeights(x))    
##})
##setGeneric("gsUniqGenes", function(obj) standardGeneric("gsUniqGenes"))
##setMethod("gsUniqGenes", "gslist", function(obj) {
##    genes <- gsGenes(obj)
##    uniqGenes <- unique(unlist(genes))
##    return(uniqGenes)
##})
##setGeneric("gsitem", function(name, desc, genes, weights) standardGeneric("gsitem"))
##setMethod("gsitem", c("character", "character", "character", "numeric"),
##          function(name, desc, genes, weights)  {
##              gs <- new("gsitem")
##              gsName(gs) <- name
##              gsDesc(gs) <- desc
##              gsGenes(gs) <- genes
##              gsWeights(gs) <- weights
##              return(gs)
##          })
##setMethod("gsitem", c("character", "character", "character", "missing"),
##          function(name, desc, genes)  {
##              return(gsitem(name, desc, genes, rep(1L, length(genes))))
##          })
##
##setMethod("show", "gsitem", function(object) {
##    genes <- gsGenes(object); glen <- length(genes); wt <- gsHasWeights(object)
##    str <- sprintf("%sGeneSet %s [Description: %s]\n%d Gene%s: %s\n",
##                   ifelse(wt, "Weighted ", ""),
##                   gsName(object), gsDesc(object),
##                   glen,
##                   ifelse(glen>1, "s", ""),
##                   ifelse(glen>5,
##                          paste(paste(genes[1:2],collapse=","),
##                                ",...,", genes[glen],sep=""),
##                          paste(genes, collapse=",")))
##    cat(str)
##})
##setMethod("show", "gslist", function(object) {
##    ol <- length(object)
##    cat("A List of", ol, "Gene Sets\n")
##    if(ol<=3) {
##        sapply(object, show)
##    } else {
##        sapply(object[1:2], show)
##        cat("...\n")
##        show(object[[ol]])
##    }
##})
##

readGmt <- function(filename) {
  stopifnot(file.exists(filename))
  lines <- readLines(filename)
  splitLines <- strsplit(lines, "\t")
  isValid <- sapply(splitLines, function(x) length(x)>=3)
  validLines <- splitLines[isValid]
  res <- lapply(validLines, function(x)  {
                    list(name=x[1], desc=x[2], genes=x[3:length(x)])
                })
  names(res) <- sapply(res, function(x) x$name)
  ## res <- .Call("read_gmt", filename, PACKAGE="BioQC")
  class(res) <- "gmtlist"
  return(res)
}

#' Convert gmtlist into a list of signed genesets
#'
#' @param gmtlist A gmtlist object, probably read-in by \code{readGmt}
#' @param posPattern Regular expression pattern of positive gene sets. It is trimmed from the original name to get the stem name of the gene set. See examples below.
#' @param negPattern Regular expression pattern of negative gene sets. It is trimmed from the original name to get the stem name of the gene set. See examples below.
#'
#' @return A list of signed genesets
#'
#' Gene set names are detected whether they are positive or negative. If neither positive nor negative 
#' @example
#' testInputList <- list(list(name="GeneSetA_UP",genes=LETTERS[1:3]),
#'                list(name="GeneSetA_DN", genes=LETTERS[4:6]),
#'                list(name="GeneSetB", genes=LETTERS[2:4]),
#'                list(name="GeneSetC_DN", genes=LETTERS[1:3]),
#'                list(name="GeneSetD_UP", genes=LETTERS[1:3]))
#' testOutputList.ignore <- gmtlist2signedGenesets(testInputList, nomatch="ignore")
#' testOutputList.pos <- gmtlist2signedGenesets(testInputList, nomatch="pos")
#' testOutputList.neg <- gmtlist2signedGenesets(testInputList, nomatch="neg")
gmtlist2signedGenesets <- function(gmtlist, posPattern="_UP$", negPattern="_DN$",
                                   nomatch=c("ignore", "pos", "neg")) {
    nomatch <- match.arg(nomatch)
    genes <- lapply(gmtlist, function(x) x$genes)
    names <- sapply(gmtlist, function(x) x$name)
    isPos <- grepl(posPattern, names)
    isNeg <- grepl(negPattern, names)
    stemNames <- gsub(negPattern, "", gsub(posPattern, "", names))
    stemFactor <- factor(stemNames, levels=unique(stemNames))
    res <- tapply(seq(along=gmtlist), stemFactor, function(i) {
                      pos <- NULL
                      neg <- NULL
                      if(!any(isPos[i]) && !any(isNeg[i])) {
                          otherGenes <- unique(unlist(genes[i]))
                          if (nomatch=="pos") {
                              pos <- otherGenes
                          } else if (nomatch=="neg") {
                              neg <- otherGenes
                          } else if (nomatch!="ignore") {
                              stop("should not be here")
                          }
                      } else {
                          if(any(isPos[i])) {
                              pos <- unique(unlist(genes[i][isPos[i]]))
                          }
                          if(any(isNeg[i])) {
                              neg <- unique(unlist(genes[i][isNeg[i]]))
                          }
                      }
                      return(list(pos=pos, neg=neg))
                  })
    class(res) <- "signed_genesets"
    return(res)
}
