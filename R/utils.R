absLog10p <- function(x) abs(log10(x))

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
