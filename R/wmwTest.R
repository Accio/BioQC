MCCORE <- 6L

wmw.test <- function(x, sub, alternative=c("two.sided", "less", "greater"), statistic=FALSE) {
  if(!any(sub)) return(ifelse(statistic, 0, 1))
  wt <- wilcox.test(x[sub],
                    x[!sub], alternative=alternative)
  return(ifelse(statistic, wt$statistic, wt$p.value))
}

setGeneric("wmwTest",function(object, sub, alternative, statistic) standardGeneric("wmwTest"))
setMethod("wmwTest", signature=c("numeric", "logical", "character", "logical"),
          function(object, sub, alternative,statistic) {
            wmw.test(object, sub, alternative=alternative, statistic=statistic)
          })
setMethod("wmwTest", signature=c("numeric", "numeric", "character", "logical"),
          function(object, sub, alternative,statistic) {
            sel <- rep(FALSE, length(object))
            sel[sub] <- TRUE
            wmw.test(object, sel, alternative=alternative, statistic=statistic)
          })
setMethod("wmwTest", signature=c("matrix", "logical", "character", "logical"),
          function(object, sub, alternative,statistic) {
            apply(object, 2L, wmw.test, sub=sub, alternative=alternative, statistic=statistic)
          })
setMethod("wmwTest", signature=c("matrix", "numeric", "character", "logical"),
          function(object, sub, alternative, statistic) {
            sel <- rep(FALSE, nrow(object))
            sel[sub] <- TRUE
            apply(object, 2L, wmw.test, sub=sel, alternative=alternative, statistic=statistic)
          })

setMethod("wmwTest", signature=c("ExpressionSet", "logical", "character", "logical"),
          function(object, sub, alternative,statistic) {
            apply(exprs(object), 2L, wmw.test, sub=sub, alternative=alternative, statistic=statistic)
          })
setMethod("wmwTest", signature=c("ExpressionSet", "numeric", "character", "logical"),
          function(object, sub, alternative, statistic) {
            sel <- rep(FALSE, nrow(object))
            sel[sub] <- TRUE
            apply(exprs(object), 2L, wmw.test, sub=sel, alternative=alternative, statistic=statistic)
          })

setMethod("wmwTest", signature=c("ExpressionSet", "list", "character", "logical"),
          function(object, sub, alternative, statistic) {
            isGmtList <- all(sapply(sub, function(x) all(c("genes", "name") %in% names(x))))
            if(!isGmtList) {
              res.raw <- lapply(sub, function(x) wmwTest(object=object, sub=x, alternative=alternative, statistic=statistic))
              res.names <- names(sub)
            } else {
              if(!"GeneSymbol" %in% colnames(fData(object))) stop("ExpressionSet must has 'GeneSymbol' as fData column\n")
              gb <- as.character(fData(object)[,"GeneSymbol"])
              oexp <- exprs(object)
              myfun <- function(x) wmwTest(oexp, gb %in% x$genes, alternative=alternative, statistic=statistic)
              
              if(require(multicore)) {
                res.raw <- mclapply(sub, myfun, mc.cores=MCCORE)
              } else {
                res.raw <- lapply(sub,myfun)
              }
              res.names <- sapply(sub, function(x) x$name)
            }
            res <- data.matrix(do.call(rbind, res.raw))
            rownames(res) <- res.names
            colnames(res) <- sampleNames(object)
            return(res)
          })

## short-cuts for default values
setMethod("wmwTest", signature=c("ANY", "ANY", "ANY", "missing"),
          function(object, sub, alternative, statistic) {
            wmwTest(object, sub, alternative, statistic=FALSE)
          })
setMethod("wmwTest", signature=c("ANY", "ANY", "missing", "ANY"),
          function(object, sub, alternative, statistic) {
            wmwTest(object, sub, alternative="two.sided", statistic)
          })
setMethod("wmwTest", signature=c("ANY", "ANY", "missing", "missing"),
          function(object, sub, alternative, statistic) {
            wmwTest(object, sub, alternative="two.sided", statistic=FALSE)
          })
