gini <- function(x)  {
    storage.mode(x) <- "double"
    isVec <- !is.matrix(x)
    if(isVec) {
        res <- .Call("gini_numeric", x, length(x))
    } else {
        res <- .Call("gini_matrix", x, nrow(x), ncol(x))
    }
    return(res)
}
