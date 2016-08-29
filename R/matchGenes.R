
##----------------------------------------##
## matchGenes for GmtList
##----------------------------------------##
matchGeneExclNA <- function(inputGene, allGene) {
    if(is.null(inputGene)) return(NULL)
    res <- match(inputGene, allGene)
    res[is.na(inputGene)] <- NA ## make sure that NA does not match to NA
    return(res)
}
#' Match genes in a GmtList to a vector of genesymbols
#'
#' @param gmtList A GmtList object
#' @param geneSymbols Gene symbols to be matched; they can come from a column in an eSet object, for example
matchGenes.default <- function(gmtList,geneSymbols) {
    if(!is(gmtList, "GmtList"))
        stop(paste("gmtlist be must of class GmtList; now it is", class(gmtList)))
    if(!is.character(geneSymbols))
        stop(paste("geneSymbols be must characters; now it is", class(geneSymbols)))
    genes <- lapply(gmtList, function(x) x$genes)
    names(genes) <- sapply(gmtList, function(x) x$name)
    indList <- sapply(genes, matchGeneExclNA, allGene=geneSymbols)
    res <- IndexList(indList)
    return(res)
}

setMethod("matchGenes", c("GmtList", "character"), function(list, object) {
              matchGenes.default(list, object)
          })
setMethod("matchGenes", c("GmtList", "matrix"), function(list, object) {
              if(is.null(rownames(object)))
                  stop("When used to map genes in GmtList directly to rows in matrix, the matrix's row names must be gene symbols")
              symbols <- rownames(object)
              matchGenes.default(list, as.character(symbols))
          })
setMethod("matchGenes", c("GmtList", "eSet"), function(list, object, col="GeneSymbol") {
              if(!is.null(col) && !col %in% colnames(fData(object)))
                  stop("When used to map genes in GmtList directly to rows in an eSet, col must be either NULL (mapped to feature names) or a column in the fData(eset)")
              if(is.null(col)) {
                  symbols <- featureNames(object)
              } else {
                  symbols <- fData(object)[,col]
              }
              matchGenes.default(list, as.character(symbols))
          })


##----------------------------------------##
## matchGenes for SignedGenesets
##----------------------------------------##
#' Match genes in a SignedGenesets to a vector of genesymbols
#'
#' @param gmtList A SignedGenesets object
#' @param geneSymbols Gene symbols to be matched; they can come from a column in an eSet object, for example

matchGenes.signedDefault <- function(signedGenesets, geneSymbols) {
    if(!is(signedGenesets, "SignedGenesets"))
        stop(paste("signedGenesets be must of class SignedGenesets; now it is", class(signedGenesets)))
    resList <- lapply(signedGenesets, function(geneset) {
        pos <- geneset$pos
        neg <- geneset$neg
        posInd <- matchGeneExclNA(geneset$pos, geneSymbols)
        negInd <- matchGeneExclNA(geneset$neg, geneSymbols)
        return(list(pos=posInd, neg=negInd))
    })
    names(resList) <- sapply(signedGenesets, function(x) x$name)
    res <- SignedIndexList(resList)
    return(res)
}

setMethod("matchGenes", c("SignedGenesets", "character"), function(list, object) {
              matchGenes.signedDefault(list, object)
          })
setMethod("matchGenes", c("SignedGenesets", "matrix"), function(list, object) {
              if(is.null(rownames(object)))
                  stop("When used to map genes in GmtList directly to rows in matrix, the matrix's row names must be gene symbols")
              symbols <- rownames(object)
              matchGenes.signedDefault(list, as.character(symbols))
          })
setMethod("matchGenes", c("SignedGenesets", "eSet"), function(list, object, col="GeneSymbol") {
              if(!is.null(col) && !col %in% colnames(fData(object)))
                  stop("When used to map genes in GmtList directly to rows in an eSet, col must be either NULL (mapped to feature names) or a column in the fData(eset)")
              if(is.null(col)) {
                  symbols <- featureNames(object)
              } else {
                  symbols <- fData(object)[,col]
              }
              matchGenes.signedDefault(list, as.character(symbols))
          })
