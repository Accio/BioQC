gini <- function(x, na.rm=FALSE)  {
  if(na.rm) x <- x[!is.na(x)]
  .Call("gini", as.numeric(x))
}
