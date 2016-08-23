wmw.test <- function(x, sub, alternative=c("two.sided", "less", "greater"), statistic=FALSE) {
  if(!any(sub)) return(ifelse(statistic, 0, 1))
  wt <- wilcox.test(x[sub],
                    x[!sub], alternative=alternative, exact=FALSE)
  return(ifelse(statistic, wt$statistic, wt$p.value))
}

## type2int and formatMatrixInd are helper functions for wmwTest and wmwSignedTest
TYPE_CODES <- c("greater"=0L, "less"=1L,
                "two.sided"=2L, "U"=3L,
                "abs.log10.greater"=4L,
                "log10.less"=5L,
                "abs.log10.two.sided"=6L,
                "Q"=7L)

type2int <- function(type) {
    if(!type %in% names(TYPE_CODES))
        stop("Should not happen! This is the wrong time")
    return(TYPE_CODES[type])
}

formatMatrixInd <- function(x, ind.list) {
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
                           if(any(x>nrow(matrix)))
                               stop("Indices out of range: they cannot exceed the row of the matrix!")
                           return(x-1L)
                       } else {
                           stop("index must be either integer vector")
                       }
                   })

    res <- list(matrix=matrix,
                indC=indC,
                isMatVec=isMatVec,
                isIndVec=isIndVec,
                rownames=names(ind.list),
                colnames=colnames(matrix))
    return(res)
  
}
wmwTest <- function(x, ind.list,
                    alternative=c("greater", "less", "two.sided", "U",
                      "abs.log10.greater","log10.less","abs.log10.two.sided","Q"), simplify=TRUE) {

    input <- formatMatrixInd(x, ind.list)

    matrix <- input$matrix
    indC <- input$indC
    typeInt <- type2int(match.arg(alternative))

    res <- .Call("wmw_test", matrix, indC, typeInt)
    
    rownames(res) <- input$rownames
    colnames(res) <- input$colnames
    
    if(simplify) {
        if(input$isMatVec & input$isIndVec) {
            res <- res[1L, 1L] 
        } else if (input$isMatVec) {
            res <- res[,1L]
        } else if (input$isIndVec) {
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
