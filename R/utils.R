absLog10p <- function(x) abs(log10(x))

filterPmat <- function(x, threshold) {
  if(missing(threshold) || is.null(threshold) || is.na(threshold) || threshold==0)
    return(x)
  fil <- apply(x, 1L, function(x) any(x <= threshold))
  x[fil,,drop=FALSE]
}
