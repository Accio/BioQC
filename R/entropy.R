atomEntropy <- function(x) ifelse(x==0, 0, x*log2(x))
atomDiversity <- function(x) {
  rel <- x/sum(x, na.rm=TRUE)
  -sum(atomEntropy(rel))
}

#' Shannon entropy
#' 
#' @param vector A vector of numbers, or characters. Discrete probability of each item is calculated and the Shannon entropy is returned.
#' @return Shannon entropy
#' 
#' Shannon entropy can be used as measures of gene expression
#' specificity, as well as measures of tissue diversity and
#' specialization. See references below.
#' 
#' We use \code{2} as base for the entropy calculation, because in this
#' base the unit of entropy is \emph{bit}.
#'
#' @references   Martinez and Reyes-Valdes (2008) Defining diversity, 
#' specialization, and gene specificity in transcriptomes through information 
#' theory. PNAS 105(28):9709--9714
#'  
#' @author Jitao David Zhang <jitao_david.zhang@roche.com>
#' 
#' @examples 
#' myVec0 <- 1:9
#' entropy(myVec0) ## log2(9)
#' myVec1 <- rep(1, 9)
#' entropy(myVec1)
#' 
#' entropy(LETTERS)
#' entropy(rep(LETTERS, 5))
#' @export
entropy <- function(vector) {
  -sum(atomEntropy(table(vector)/length(vector)))
}

#' Entropy-based gene-expression specificity
#' @param mat A matrix (usually an expression matrix), with genes (features) in rows and samples in columns.
#' @param norm Logical, whether the specificity should be normalized by \code{log2(ncol(mat))}.
#' @return A vector of the length of the row number of the input matrix, namely the specificity score of genes.
#' 
#' @references   Martinez and Reyes-Valdes (2008) Defining diversity, 
#' specialization, and gene specificity in transcriptomes through information 
#' theory. PNAS 105(28):9709--9714
#' 
#' @seealso \code{\link{entropy}}
#' @examples
#' myMat <- rbind(c(3,4,5),c(6,6,6), c(0,2,4))
#' entropySpecificity(myMat)
#' entropySpecificity(myMat, norm=TRUE)
#'
#' myRandomMat <- matrix(runif(1000), ncol=20)
#' entropySpecificity(myRandomMat)
#' entropySpecificity(myRandomMat, norm=TRUE)
#' @export
entropySpecificity <- function(mat, norm=FALSE) {
  mat.rel <- apply(mat, 2L, function(x) x/sum(x, na.rm=TRUE))
  speci <- apply(mat.rel, 1L, function(x) {
    gm <- mean(x, na.rm=TRUE)
    rel <- x/gm
    mean(atomEntropy(rel), na.rm=TRUE)
  })
  if(norm)
    speci <- speci/log2(ncol(mat))
  return(speci)
}

#' Entropy-based sample diversity
#' @param mat A matrix (usually an expression matrix), with genes (features) in rows and samples in columns.
#' @param norm Logical, whether the diversity should be normalized by \code{log2(nrow(mat))}.
#' @return A vector as long as the column number of the input matrix
#' 
#' @references   Martinez and Reyes-Valdes (2008) Defining diversity, 
#' specialization, and gene specificity in transcriptomes through information 
#' theory. PNAS 105(28):9709--9714
#' 
#' @seealso \code{\link{entropy}} and \code{\link{sampleSpecialization}}
#' @examples
#' myMat <- rbind(c(3,4,5),c(6,6,6), c(0,2,4))
#' entropyDiversity(myMat)
#' entropyDiversity(myMat, norm=TRUE)
#'
#' myRandomMat <- matrix(runif(1000), ncol=20)
#' entropyDiversity(myRandomMat)
#' entropyDiversity(myRandomMat, norm=TRUE)
#' @export
entropyDiversity <- function(mat, norm=FALSE) {
  res <- apply(mat, 2L, atomDiversity)
  if(norm)
    res <- res/log2(nrow(mat))
  res
}


#' Entropy-based sample specialization
#' @param mat A matrix (usually an expression matrix), with genes (features) in rows and samples in columns.
#' @param norm Logical, whether the specialization should be normalized by \code{log2(ncol(mat))}.
#' @return A vector as long as the column number of the input matrix
#' 
#' @references   Martinez and Reyes-Valdes (2008) Defining diversity, 
#' specialization, and gene specificity in transcriptomes through information 
#' theory. PNAS 105(28):9709--9714
#' 
#' @seealso \code{\link{entropy}} and \code{\link{entropyDiversity}}
#' @examples
#' myMat <- rbind(c(3,4,5),c(6,6,6), c(0,2,4))
#' sampleSpecialization(myMat)
#' sampleSpecialization(myMat, norm=TRUE)
#'
#' myRandomMat <- matrix(runif(1000), ncol=20)
#' sampleSpecialization(myRandomMat)
#' sampleSpecialization(myRandomMat, norm=TRUE)
#' @export
sampleSpecialization <- function(mat, norm=TRUE) {
  speci <- entropySpecificity(mat, norm=FALSE)
  rel <- apply(mat, 2L, function(x) x/sum(x, na.rm=TRUE))
  res <- apply(rel, 2L, function(x) sum(x*speci, na.rm=TRUE))
  if(norm)
    res <- res/log2(ncol(mat))
  res
}
