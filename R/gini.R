#'Calculate Gini Index of a numeric vector
#'
#'Calculate the Gini index of a numeric vector.
#'
#'@param x A numeric vector.
#'
#'@details The Gini index (Gini coefficient) is a measure of statistical dispersion. A Gini
#'coefficient of zero expresses perfect equality where all values are
#'the same. A Gini coefficient of one expresses maximal inequality among values.
#'
#'@return A numeric value between 0 and 1.
#'
#'@references Gini. C. (1912) \emph{Variability and Mutability}, C. Cuppini, Bologna
#'156 pages.
#'
#'@author Jitao David Zhang
#'
#'@examples 
#'testValues <- runif(100)
#'gini(testValues)
#'    
gini <- function(x)  {
    storage.mode(x) <- "double"
    isVec <- !is.matrix(x)
    if(isVec) {
        x <- x[!is.na(x)]
        x <- sort(x, decreasing=FALSE)
        res <- .Call(C_gini_numeric, x, length(x))
    } else {
        res <- .Call(C_gini_matrix, x, nrow(x), ncol(x))
    }
    return(res)
}
