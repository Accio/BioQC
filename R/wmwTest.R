wmw.test <- function(x, sub, alternative=c("two.sided", "less", "greater"), statistic=FALSE) {
  if(!any(sub)) return(ifelse(statistic, 0, 1))
  wt <- wilcox.test(x[sub],
                    x[!sub], alternative=alternative, exact=FALSE)
  return(ifelse(statistic, wt$statistic, wt$p.value))
}

wmwTest <- function(x, ind.list,
                    alternative=c("greater", "less", "two.sided", "U",
                      "abs.log10.greater","log10.less","abs.log10.two.sided","Q"), simplify=TRUE, useNonFactor=FALSE) {

  isMatVec <- FALSE
  isIndVec <- FALSE
  if(is(x, "eSet")) {
    matrix <- exprs(x)
  } else if (!is.matrix(x) & is.numeric(x)) {
    matrix <- matrix(x, ncol=1L)
    isMatVec <- TRUE
  } else if (is.matrix(x)) {
    matrix <- x
  } else {
    stop("'matrix' must be a numeric matrix, or a numeric vector, or an eSet object")
  }

  if(storage.mode(matrix)!="double")
      storage.mode(matrix) <- "double"
  
  if(is.numeric(ind.list) || is.logical(ind.list)) {
    ind.list <- list(ind.list)
    isIndVec <- TRUE
  } else if (is(ind.list, "gmtlist")) {
    if(!(is(x, "eSet") & "GeneSymbol" %in% colnames(fData(x))))
      stop("When ind.list is a 'gmtlist', 'x' must be an eSet object with a column 'GeneSymbol' in fData.")
    genes <- lapply(ind.list, function(x) x$genes)
    names(genes) <- sapply(ind.list, function(x) x$name)
    ind.list <- sapply(genes, function(g) match(g, fData(x)$GeneSymbol))
  }
  
  indC <- lapply(ind.list, function(x) {
    if(is.logical(x))
      x <- which(x)
    if(is.numeric(x)) {
      x <- as.integer(x[!is.na(x)])
      if(any(x<1L))
          stop("Indices in ind.list must be equal to or greater than one!")
      return(x-1L)
    } else {
      stop("index must be either integer vector")
    }
  })
  
  type <- match.arg(alternative)
  if(type=="greater") {
    val <- 0L
  } else if (type=="less") {
    val <- 1L
  } else if (type=="two.sided") {
    val <- 2L
  } else if (type=="U") {
    val <- 3L
  } else if (type=="abs.log10.greater") {
    val <- 4L
  } else if (type=="log10.less") {
    val <- 5L
  } else if (type=="abs.log10.two.sided") {
    val <- 6L
  } else if (type=="Q") {
    val <- 7L
  } else {
    stop("Should not happen")
  }
  if(useNonFactor) {
    res <- .Call("wmw_test", indC, matrix, val)
  } else {
    res <- .Call("wmw_test_nonfactor", indC, matrix, val)
  } 
  rownames(res) <- names(ind.list)
  colnames(res) <- colnames(matrix)
  if(simplify) {
    if(isMatVec & isIndVec) {
      res <- res[1L, 1L] 
    } else if (isMatVec) {
      res <- res[,1L]
    } else if (isIndVec) {
      res <- res[1L,]
    }
  }
  
  return(res)
}

##setGeneric("wmwTest",function(object, sub, alternative, statistic) standardGeneric("wmwTest"))
##setGeneric("wmwTest",function(exprs, index, alternative) standardGeneric("wmwTest"))
##setMethod("wmwTest", signature=c("ANY", "ANY", "character") , function(exprs, index, alternative) {
##  wmwTestC(exprs, index, alternative=alternative)
##})
##setMethod("wmwTest", signature=c("eSet", "gmtlist", "character"), function(exprs, index, alternative) {
##  if(!"GeneSymbol" %in% colnames(fData(exprs)))
##    stop("ExpressionSet must has 'GeneSymbol' as fData column which contains gene symbols used in the GMT files\n")
##  gb <- as.character(fData(exprs)[,"GeneSymbol"])
##  ind <- lapply(index, function(x) {
##    rind <- match(x$genes, gb)
##    rind <- rind[!is.na(rind)]
##    return(rind-1L)
##  })
##  names(ind) <- sapply(index, function(x) x$name)
##  wmwTest(exprs(exprs), ind, alternative=alternative)
##})
