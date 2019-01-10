
##----------------------------------------##
## matchGenes for GmtList
##----------------------------------------##
matchGeneExclNA <- function(inputGene, allGene) {
    if(is.null(inputGene)) return(NULL)
    res <- match(inputGene, allGene)
    res[is.na(inputGene)] <- NA ## make sure that NA does not match to NA
    return(res)
}

matchGenes.default <- function(gmtList,geneSymbols) {
    if(!is(gmtList, "GmtList"))
        stop(paste("gmtlist be must of class GmtList; now it is", class(gmtList)))
    if(!is.character(geneSymbols))
        stop(paste("geneSymbols be must characters; now it is", class(geneSymbols)))
    genes <- lapply(gmtList, function(x) x$genes)
    names(genes) <- sapply(gmtList, function(x) x$name)
    indList <- lapply(genes, matchGeneExclNA, allGene=geneSymbols)
    res <- IndexList(indList)
    return(res)
}
#' Match genes in a list-like object to a vector of genesymbols
#'
#' @param list A GmtList, list, character or SignedGenesets object
#' @param object Gene symbols to be matched; they can come from a vector of character strings, or
#' a column in the fData of an \code{eSet} object.
#' @param ... additional arguments like \code{col}
#' @param col Column name of \code{fData} in an \code{eSet} object, or \code{genes} in an \code{DGEList} object, to specify where gene symbols are stored.
#' The default value is set to "GeneSymbol"
#' @name matchGenes
#' @examples
#' ## test GmtList, character
#' testGenes <- sprintf("gene%d", 1:10)
#' testGeneSets <- GmtList(list(gs1=c("gene1", "gene2"), gs2=c("gene9", "gene10"), gs3=c("gene100")))
#' matchGenes(testGeneSets, testGenes)
#' 
#' ## test GmtList, matrix
#' testGenes <- sprintf("gene%d", 1:10)
#' testGeneSets <- GmtList(list(gs1=c("gene1", "gene2"), gs2=c("gene9", "gene10"), gs3=c("gene100")))
#' testGeneExprs <- matrix(rnorm(100), nrow=10, dimnames=list(testGenes, sprintf("sample%d", 1:10)))
#' matchGenes(testGeneSets, testGeneExprs)
#' 
#' ## test GmtList, eSet
#' testGenes <- sprintf("gene%d", 1:10)
#' testGeneSets <- GmtList(list(gs1=c("gene1", "gene2"), gs2=c("gene9", "gene10"), gs3=c("gene100")))
#' testGeneExprs <- matrix(rnorm(100), nrow=10, dimnames=list(testGenes, sprintf("sample%d", 1:10)))
#' testFeat <- data.frame(GeneSymbol=rownames(testGeneExprs), row.names=testGenes)
#' testPheno <- data.frame(SampleId=colnames(testGeneExprs), row.names=colnames(testGeneExprs))
#' testEset <- ExpressionSet(assayData=testGeneExprs,
#'     featureData=AnnotatedDataFrame(testFeat),
#'     phenoData=AnnotatedDataFrame(testPheno))
#' matchGenes(testGeneSets, testGeneExprs)
#' ## force using row names
#' matchGenes(testGeneSets, testEset, col=NULL)
#' 
#'  ## test GmtList, DGEList
#'  if(requireNamespace("edgeR")) {
#'     mat <- matrix(rnbinom(10000, mu=5, size=2), ncol=4)
#'     rownames(mat) <- sprintf("gene%d", 1:nrow(mat))
#'     y <- edgeR::DGEList(counts=mat, group=rep(1:2, each=2))
#'
#'     ## if genes are not set, row names of the count matrix will be used for lookup
#'     myGeneSet <- GmtList(list(gs1=rownames(mat)[1:3], gs2=rownames(mat)[5:10]))
#'     matchGenes(myGeneSet, y)
#'
#'     ## alternatively, use 'col' parameter to specify the column in 'genes'
#'     y2 <- edgeR::DGEList(counts=mat,
#'       group=rep(1:2, each=2),
#'       genes=data.frame(GeneIdentifier=rownames(mat), row.names=rownames(mat)))
#'     matchGenes(myGeneSet, y2, col="GeneIdentifier")
#'  }
#' 
#' ## test character, character
#' matchGenes(c("gene1", "gene2"), testGenes)
#' 
#' ## test character, matrix
#' matchGenes(c("gene1", "gene2"), testGeneExprs)
#' 
#' ## test character, eset
#' matchGenes(c("gene1", "gene2"), testEset)
NULL

#'@rdname matchGenes
setMethod("matchGenes", c("GmtList", "character"), function(list, object) {
              matchGenes.default(list, object)
          })

#'@rdname matchGenes
setMethod("matchGenes", c("GmtList", "matrix"), function(list, object) {
              if(is.null(rownames(object)))
                  stop("When used to map genes in GmtList directly to rows in matrix, the matrix's row names must be gene symbols")
              symbols <- rownames(object)
              matchGenes.default(list, as.character(symbols))
          })


#'@rdname matchGenes
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

#'@rdname matchGenes
setMethod("matchGenes", c("character", "character"), function(list, object) {
              tempList <- GmtList(list(TempGeneSet=list))
              matchGenes.default(tempList, object)
          })
#'@rdname matchGenes
setMethod("matchGenes", c("character", "matrix"), function(list, object) {
              tempList <- GmtList(list(TempGeneSet=list))
              matchGenes(tempList, object)
          })
#'@rdname matchGenes
setMethod("matchGenes", c("character", "eSet"), function(list, object) {
              tempList <- GmtList(list(TempGeneSet=list))
              matchGenes(tempList, object)
          })

#'@rdname matchGenes
setMethod("matchGenes", c("GmtList", "DGEList"), function(list, object, col="GeneSymbol") {
  if(is.null(object$genes)) {
    col <- NULL
  } else {
    if (!is.null(col) && !col %in% colnames(object$genes)) {
      stop("'col' not found in the column names of object$genes")
    }
  }
  
  if (is.null(col)) {
    symbols <- rownames(object$counts)
  } else {
    symbols <- object$genes[, col]
  }
  matchGenes.default(list, as.character(symbols))
})

##----------------------------------------##
## matchGenes for SignedGenesets
##----------------------------------------##


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
#'@rdname matchGenes
setMethod("matchGenes", c("SignedGenesets", "character"), function(list, object) {
              matchGenes.signedDefault(list, object)
          })
#'@rdname matchGenes
setMethod("matchGenes", c("SignedGenesets", "matrix"), function(list, object) {
              if(is.null(rownames(object)))
                  stop("When used to map genes in GmtList directly to rows in matrix, the matrix's row names must be gene symbols")
              symbols <- rownames(object)
              matchGenes.signedDefault(list, as.character(symbols))
          })
#'@rdname matchGenes
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
