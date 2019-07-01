#' Absolute base-10 logarithm of p-values
#' 
#' @param x Numeric vector or matrix
#' 
#' The function returns the absolute values of base-10 logarithm of p-values.
#' 
#' @details 
#' The logarithm transformation of p-values is commonly used to visualize
#' results from statistical tests. Although it may cause misunderstanding
#' and therefore its use is disapproved by some experts, it helps to
#' visualize and interpret results of statistical tests intuitively.
#' 
#' The function transforms p-values with base-10 logarithm, and returns its
#' absolute value. The choice of base 10 is driven by the simplicity of
#' interpreting the results.
#' 
#' @return Numeric vector or matrix.
#' @author Jitao David Zhang <jitao_david.zhang@roche.com>
#' @importFrom utils head
#' @examples 
#' testp <- runif(1000, 0, 1)
#' testp.al <- absLog10p(testp)
#' 
#' print(head(testp))
#' print(head(testp.al))
#' 
#' @export
absLog10p <- function(x) abs(log10(x))

#' Filter rows of p-value matrix under the significance threshold
#' 
#' @param x A matrix of p-values. It must be raw p-values and should not be transformed (e.g. logarithmic).
#' @param threshold A numeric value, the minimal p-value used to filter
#' rows. If missing, given the values of \code{NA}, \code{NULL} or
#' number \code{0}, no filtering will be done and the input matrix will
#' be returned.
#' 
#' @return   Matrix of p-values. If no line is left, a empty matrix of the same dimension as input will be returned.
#' @author Jitao David Zhang <jitao_david.zhang@roche.com>
#' 
#' @examples 
#' set.seed(1235)
#' testMatrix <- matrix(runif(100,0,1), nrow=10)
#' 
#' ## filtering
#' (testMatrix.filter <- filterPmat(testMatrix, threshold=0.05))
#' ## more strict filtering
#' (testMatrix.strictfilter <- filterPmat(testMatrix, threshold=0.01))
#' ## no filtering
#' (testMatrix.nofilter <- filterPmat(testMatrix))
#' @export 
filterPmat <- function(x, threshold) {
  if(missing(threshold) || is.null(threshold) || is.na(threshold) || threshold==0)
    return(x)
  fil <- apply(x, 1L, function(x) any(x <= threshold))
  x[fil,,drop=FALSE]
}

#' Simplify matrix in case of single row/columns
#'
#' @param matrix A matrix of any dimension
#'
#' If only one row/column is present, the dimension is dropped and a vector will be returned
#'
#' @examples
#' testMatrix <- matrix(round(rnorm(9),2), nrow=3)
#' simplifyMatrix(testMatrix)
#' simplifyMatrix(testMatrix[1L,,drop=FALSE])
#' simplifyMatrix(testMatrix[,1L,drop=FALSE])
#' @export
simplifyMatrix <- function(matrix) {
    if(nrow(matrix)==1L & ncol(matrix)==1L)  {
        matrix <- matrix[1L,1L]
    } else if (nrow(matrix)==1L) {
        matrix <- matrix[1L,]
    } else if (ncol(matrix)==1L)  {
        matrix <- matrix[,1L]
    } 
    return(matrix)
}
