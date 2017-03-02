gini <- function(x)  {
    isVec <- !is.matrix(x)
    if(isVec)
        x <- matrix(x, nrow=1L, byrow=TRUE)
    storage.mode(x) <- "double"
    res <- .Call("gini_matrix", x, nrow(x), ncol(x))
    if(isVec)
        res <- res[1]
    return(res)
}
