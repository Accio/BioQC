#' Getting leading-edge indices from a vector
#' @param x A numeric vector (\code{getLeadingEdgeIndexFromVector}) or a numeric matrix (\code{getLeadingEdgeIndexFromMatrix}).
#' @param index An integer vector, indicating the indices of genes in a gene-set.
#' @param comparison Character string, are values greater than or less than the reference value considered as leading-edge? This depends on the type of value requested by the user in \code{wmwTest}.
#' @param reference Character string, which reference is used? If \code{background}, genes with expression higher than the median of the background are reported. Otherwise in the case of \code{geneset}, genes with expression higher than the median of the gene-set is reported. Default is \code{background}, which is consistent with the results of the Wilcoxon-Mann-Whitney tests.
#' 
#' @return An integer vector, indicating the indices of leading-edge genes.
#' @seealso \code{\link{wmwTest}}
#' @examples
#' myProfile <- c(rnorm(5, 3), rnorm(15, -3), rnorm(100, 0))
#' getLeadingEdgeIndexFromVector(myProfile, 1:20)
#' getLeadingEdgeIndexFromVector(myProfile, 1:20, comparison="less")
#' getLeadingEdgeIndexFromVector(myProfile, 1:20, comparison="less", reference="geneset")
#' myProfile2 <- c(rnorm(15, 3), rnorm(5, -3), rnorm(100, 0))
#' myProfileMat <- cbind(myProfile, myProfile2)
#' getLeadingEdgeIndexFromMatrix(myProfileMat, 1:20)
#' getLeadingEdgeIndexFromMatrix(myProfileMat, 1:20, comparison="less")
#' getLeadingEdgeIndexFromMatrix(myProfileMat, 1:20, comparison="less", reference="geneset")
#' @export
getLeadingEdgeIndexFromVector <- function(x, index,
                                          comparison = c("greater", "less"),
                                          reference = c("background", "geneset")) {
  indexValues <- x[index]
  reference <- match.arg(reference)
  comparison <- match.arg(comparison)
  ref <- switch(reference,
                geneset = median(indexValues, na.rm=TRUE),
                background = median(x, na.rm=TRUE))
  isLeading <- switch(comparison,
                      greater = indexValues >= ref,
                      less = indexValues <= ref)
  res <- index[isLeading]
  return(res)
}

#' @describeIn getLeadingEdgeIndexFromVector \code{x} is a \code{matrix}.
#' @export
getLeadingEdgeIndexFromMatrix <- function(x, index, 
                                          comparison = c("greater", "less"),
                                          reference = c("background", "geneset")) {
  reference <- match.arg(reference)
  comparison <- match.arg(comparison)
  res <- apply(x, 2, getLeadingEdgeIndexFromVector, index=index,
        reference=reference, comparison=comparison)
  if(is.matrix(res)) {
    resList <- lapply(seq_len(ncol(res)), function(i) res[,i,drop=TRUE])
    names(resList) <- colnames(res)
    res <- resList
  }
  return(res)
}


#' Identify BioQC leading-edge genes of one gene-set
#' 
#' @param matrix A numeric matrix
#' @param indexVector An integer vector, giving indices of a gene-set of interest
#' @param valType Value type, consistent with the types in \code{wmwTest}
#' @param thr Threshold of the value, greater or less than which the gene-set is considered significantly enriched in one sample
#' @param reference Character string, which reference is used? If \code{background}, genes with expression higher than the median of the background are reported. Otherwise in the case of \code{geneset}, genes with expression higher than the median of the gene-set is reported. Default is \code{background}, which is consistent with the results of the Wilcoxon-Mann-Whitney tests.
#' @return A list of integer vectors.
#' 
#' BioQC leading-edge genes are defined as those features whose expression is higher than the median expression of the background in a sample.  The function identifies leading-edge genes of a given dataset (specified by the index vector) in a number of samples (specified by the matrix, with genes/features in rows and samples in columns) in three steps. The function calls \code{wmwTest} to run BioQC and identify samples in which the gene-set is significantly enriched. The enrichment criteria is specified by \code{valType} and \code{thr}. Then the function identifies genes in the gene-set that have greater or less expresion than the median value of the \code{reference} in those samples showing significant enrichment. Finally, it reports either leading-edge genes in individual samples, or the intersection/union of leading-edge genes in multiple samples.
#' 
#' @seealso \code{\link{wmwTest}}
#' @importFrom stats median
#' @examples
#' myProfile <- c(rnorm(5, 3), rnorm(15, -3), rnorm(100, 0))
#' myProfile2 <- c(rnorm(15, 3), rnorm(5, -3), rnorm(100, 0))
#' myProfile3 <- c(rnorm(10, 5), rnorm(10, 0), rnorm(100, 0))
#' myProfileMat <- cbind(myProfile, myProfile2, myProfile3)
#' wmwLeadingEdge(myProfileMat, 1:20, valType="p.greater")
#' wmwLeadingEdge(myProfileMat, 1:20, valType="log10p.less")
#' wmwLeadingEdge(myProfileMat, 1:20, valType="U", reference="geneset")
#' wmwLeadingEdge(myProfileMat, 1:20, valType="abs.log10p.greater")
#' @export
wmwLeadingEdge <- function(matrix, 
                               indexVector, 
                               valType = c("p.greater", "p.less", "p.two.sided", "U", "abs.log10p.greater",
                                           "log10p.less", "abs.log10p.two.sided", "Q", "r", "f", "U1", "U2"),
                               thr = 0.05, 
                               reference = c("background", "geneset")) {
  if (missing(valType)) {
    valType <- "p.greater"
  } else {
    valType <- match.arg(valType)
  }
  reference <- match.arg(reference)

  lowerIsSig <- grepl("^p\\.", valType) || valType %in% c("log10p.less")
  higherIsSig <- grepl("^abs\\.", valType) || valType %in% c("U", "Q", "r", "f", "U1", "U2")
  if (!lowerIsSig & !higherIsSig) {
    stop(paste("Not implemented: is lower or higher values of", valType, "significant?"))
  }
  
  higherIsLeading <- grepl("greater", valType) || valType %in% c("U", "Q", "r", "f", "U1", "U2")
  lowerIsLeading <- grepl("less", valType)
  if (!higherIsLeading & !lowerIsLeading) {
    stop(paste("Not implemented: leading-edge genes are not implemented for ", valType, "yet"))
  } else {
    stopifnot(!(higherIsLeading & lowerIsLeading))
    compType <- ifelse(higherIsLeading, "greater", "less")
  }
  wmwTestRes <- wmwTest(matrix, indexVector, 
                                valType=valType, simplify=TRUE)
  if(lowerIsSig) {
    isSig <- wmwTestRes <= thr
  } else {
    isSig <- wmwTestRes >= thr
  }
  if(!any(isSig)) {
    warning("No single sample reported significant results. Please make sure the threshold is set correcly.")
    return(NULL)
  }
  sigMatrix <- matrix[, isSig, drop=FALSE]

  resList <- getLeadingEdgeIndexFromMatrix(sigMatrix, index=indexVector,
                                           reference=reference, 
                                           comparison=compType)
  return(resList)
}
