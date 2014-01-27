setClass("gsitem",
         representation(name="character",
                        desc="character",
                        genes="character",
                        weights="numeric",
                        indices="integer"))
setClass("gslist",  contains="list")
setAs(from="list", to="gslist", def=function(from, to) {
    new(to, from)
})

setGeneric("gsName", function(obj) standardGeneric("gsName"))
setGeneric("gsDesc", function(obj) standardGeneric("gsDesc"))
setGeneric("gsGenes", function(obj) standardGeneric("gsGenes"))
setGeneric("gsWeights", function(obj) standardGeneric("gsWeights"))
setGeneric("gsHasWeights", function(obj) standardGeneric("gsHasWeights"))
setGeneric("gsIndices", function(obj) standardGeneric("gsIndices"))
setMethod("gsName", "gsitem", function(obj) obj@name)
setMethod("gsDesc", "gsitem", function(obj) obj@desc)
setMethod("gsGenes", "gsitem", function(obj) obj@genes)
setMethod("gsWeights", "gsitem", function(obj) obj@weights)
setMethod("gsHasWeights", "gsitem", function(obj) !all(obj@weights==1L))
setMethod("gsIndices", "gsitem", function(obj) obj@indices)
setMethod("length", "gsitem", function(x) length(x@genes)) 

setMethod("gsName", "gslist", function(obj) sapply(obj, gsName))
setMethod("gsDesc", "gslist", function(obj) sapply(obj, gsDesc))
setMethod("gsGenes", "gslist", function(obj) lapply(obj, gsGenes))
setMethod("gsWeights", "gslist", function(obj) lapply(obj, gsWeights))
setMethod("gsHasWeights", "gslist", function(obj) sapply(obj, gsHasWeights))
setMethod("gsIndices", "gslist", function(obj) sapply(obj, gsIndices))
          
setGeneric("gsName<-", function(obj,value) standardGeneric("gsName<-"))
setGeneric("gsDesc<-", function(obj,value) standardGeneric("gsDesc<-"))
setGeneric("gsGenes<-", function(obj,value) standardGeneric("gsGenes<-"))
setGeneric("gsWeights<-", function(obj,value) standardGeneric("gsWeights<-"))
setGeneric("gsIndices<-", function(obj,value) standardGeneric("gsIndices<-"))
setMethod("gsName<-", c("gsitem","character"), function(obj,value) {
    obj@name <- value[1]
    return(obj)
})
setMethod("gsDesc<-", c("gsitem", "character"), function(obj,value) {
    obj@desc <- value[1]
    return(obj)
})
setMethod("gsGenes<-", c("gsitem","character"), function(obj,value){
    obj@genes <- unique(value)
    return(obj)
})
setMethod("gsWeights<-", c("gsitem","numeric"), function(obj,value) {
    if(length(value)!=length(obj))
        stop("weights must be of the same length as the genes\n")
    obj@weights <- value
    return(obj)
})
setMethod("gsIndices<-", c("gsitem","integer"), function(obj,value) {
    if(length(value)!=length(obj))
        stop("indices must be of the same length as the genes\n")
    obj@indices <- value
    return(obj)
})

setGeneric("gsitem", function(name, desc, genes, weights) standardGeneric("gsitem"))
setMethod("gsitem", c("character", "character", "character", "numeric"),
          function(name, desc, genes, weights)  {
              gs <- new("gsitem")
              gsName(gs) <- name
              gsDesc(gs) <- desc
              gsGenes(gs) <- genes
              gsWeights(gs) <- weights
              return(gs)
          })
setMethod("gsitem", c("character", "character", "character", "missing"),
          function(name, desc, genes)  {
              return(gsitem(name, desc, genes, rep(1L, length(genes))))
          })

setMethod("show", "gsitem", function(object) {
    genes <- gsGenes(object); glen <- length(genes); wt <- gsHasWeights(object)
    str <- sprintf("%sGeneSet %s [Description: %s]\n%d Gene%s: %s\n",
                   ifelse(wt, "Weighted ", ""),
                   gsName(object), gsDesc(object),
                   glen,
                   ifelse(glen>1, "s", ""),
                   ifelse(glen>3,
                          paste(paste(genes[1:2],collapse=","),
                                ",...,", genes[gölen],sep=""),
                          paste(genes, collapse=",")))
    cat(str)
})
setMethod("show", "gslist", function(object) {
    ol <- length(object)
    cat("A List of", ol, "Gene Sets\n")
    if(ol<=3) {
        sapply(object, show)
    } else {
        sapply(object[1:2], show)
        cat("...\n")
        show(object[[ol]])
    }
})

readGmt <- function(filename) {
  stopifnot(file.exists(filename))
  filename <- path.expand(filename)
  .Call("read_gmt", filename)
}

