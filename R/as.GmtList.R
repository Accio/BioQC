#' Convert a list of gene symbols into a gmtlist
#'
#' @param list A named list with character vectors of genes. Names will become names of gene sets; character vectors will become genes
#' @param description Character, description of gene-sets. The value will be expanded to the same length of the list.
#' @param uniqGenes Logical, whether redundant genes should be made unique?
#' @param namespace Character or \code{NULL}, namespace of the gene-set
#' 
#' @examples
#' testVec <- list(GeneSet1=c("AKT1", "AKT2"),
#'                GeneSet2=c("MAPK1", "MAPK3"),
#'                GeneSet3=NULL)
#' testVecGmtlist <- as.GmtList(testVec)
#' 
#' @export
as.GmtList <- function(list, description=NULL,
                       uniqGenes=TRUE,
                       namespace=NULL) {
  if(is.null(names(list)))
    stop('The input list must have non-null names')
  names <- names(list)
  if(!is.null(description)) {
    descs <- rep_len(description, length.out=length(list))
  }
  res <- lapply(seq(along=list),
                function(i) {
                  if(!is.null(description)) {
                    desc <- descs[i]
                  } else {
                    desc <- NULL
                  }
                  genes <- list[[i]]
                  if(uniqGenes) {
                    genes <- unique(genes)
                  }
                  list(name=names[i],
                       desc=desc,
                       genes=genes,
                       namespace=namespace)
                })
  names(res) <- names
  return(GmtList(res))
}