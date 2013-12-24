atomEntropy <- function(x) ifelse(x==0, 0, x*log2(x))
atomDiversity <- function(x) {
  rel <- x/sum(x, na.rm=TRUE)
  -sum(atomEntropy(rel))
}
entropy <- function(vector) {
  -sum(atomEntropy(table(vector)/length(vector)))
}
entropyDiversity <- function(mat, norm=FALSE) {
  res <- apply(mat, 2L, atomDiversity)
  if(norm)
    res <- res/log2(nrow(mat))
  res
}
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
sampleSpecialization <- function(mat, norm=TRUE) {
  speci <- entropySpecificity(mat, norm=FALSE)
  rel <- apply(mat, 2L, function(x) x/sum(x, na.rm=TRUE))
  res <- apply(rel, 2L, function(x) sum(x*speci, na.rm=TRUE))
  if(norm)
    res <- res/log2(ncol(mat))
  res
}
