wmw.test <- function(x, sub, alternative=c("two.sided", "less", "greater"), statistic=FALSE) {
  if(!any(sub)) return(ifelse(statistic, 0, 1))
  wt <- wilcox.test(x[sub],
                    x[!sub], alternative=alternative)
  return(ifelse(statistic, wt$statistic, wt$p.value))
}

wmwTestC <- function(x, sub, alternative=c("two.sided", "less", "greater"), statistic=FALSE) {
  alternative <- match.arg(alternative)
  if(!any(sub)) return(ifelse(statistic, 0, 1))
  if(statistic) {
    val <- 3L
  } else if (alternative=="less") {
    val <- 1L
  } else if (alternative=="greater") {
    val <- 0L
  } else if (alternative=="two.sided") {
    val <- 2L
  } else {
    stop("Should not happen")
  }
  return(.Call("wmw_test", which(sub)-1L, matrix(x, ncol=1L), val=val))
}

wmwTestCind <- function(matrix, ind.list, alternative=c("greater", "less", "two.sided", "U")) {
  isMatVec <- FALSE
  isIndVec <- FALSE
  if(is(matrix, "eSet")) {
    matrix <- exprs(matrix)
  } else if (!is.matrix(matrix) & is.numeric(matrix)) {
    matrix <- matrix(matrix, ncol=1L)
    isMatVec <- TRUE
  }
  if(is.numeric(ind.list) || is.logical(ind.list)) {
    ind.list <- list(ind.list)
    isIndVec <- TRUE
  }
  
  indC <- lapply(ind.list, function(x) {
    if(is.logical(x))
      x <- which(x)
    if(is.numeric(x)) {
      x <- as.integer(x[!is.na(x)])
      return(x-1L)
    } else {
      stop("index must be either integer vector")
    }
  })
  
  if(!is.matrix(matrix) || !is.numeric(matrix))
    stop("matrix must be a numeric matrix")
  type <- match.arg(alternative)
  if(type=="greater") {
    val <- 0L
  } else if (type=="less") {
    val <- 1L
  } else if (type=="two.sided") {
    val <- 2L
  } else if (type=="U") {
    val <- 3L
  } else {
    stop("Should not happen")
  }
  res <- .Call("wmw_test", indC, matrix, val)
  rownames(res) <- names(ind.list)
  colnames(res) <- colnames(matrix)
  if(isMatVec) res <- res[,1L]
  if(isIndVec) res <- res[1L,]

  return(res)
}

##setGeneric("wmwTest",function(object, sub, alternative, statistic) standardGeneric("wmwTest"))
setGeneric("wmwTest",function(exprs, index, alternative) standardGeneric("wmwTest"))
setMethod("wmwTest", signature=c("matrix", "numeric", "character") , function(exprs, index, alternative) {
  wmwTestCind(exprs, index, alternative=alternative)
})
setMethod("wmwTest", signature=c("eSet", "gmtlist", "character"), function(exprs, index, alternative) {
  if(!"GeneSymbol" %in% colnames(fData(exprs)))
    stop("ExpressionSet must has 'GeneSymbol' as fData column which contains gene symbols used in the GMT files\n")
  gb <- as.character(fData(exprs)[,"GeneSymbol"])
  ind <- lapply(index, function(x) {
    rind <- match(x$genes, gb)
    rind <- rind[!is.na(rind)]
    return(rind-1L)
  })
  names(ind) <- sapply(index, function(x) x$name)
  wmwTest(exprs(exprs), ind, alternative=alternative)
})
