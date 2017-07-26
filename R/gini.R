gini <- function(x)  {
    storage.mode(x) <- "double"
    isVec <- !is.matrix(x)
    if(isVec) {
        x <- x[!is.na(x)]
        x <- sort(x, decreasing=FALSE)
        res <- .Call("gini_numeric", x, length(x))
    } else {
        res <- .Call("gini_matrix", x, nrow(x), ncol(x))
    }
    return(res)
}
