gini <- function(x)  {
  .Call("gini", as.numeric(x))
}
