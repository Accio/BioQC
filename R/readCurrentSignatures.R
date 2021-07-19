#' Load current BioQC signatures
#' @param uniqGenes Logical, whether duplicated genes should be removed, passed
#' to \code{\link{readGmt}}.
#' @param namespace Character, namespace of the gene-set, or code{NULL}, passed
#' to \code{\link{readGmt}}
#' @return A GmtList
#' @seealso \code{\link{readGmt}}
#' @examples
#' readCurrentSignatures()
#' @export
readCurrentSignatures <- function(uniqGenes=TRUE, namespace=NULL) {
    file <- system.file("extdata/exp.tissuemark.bioqc.roche.symbols.gmt", package="BioQC")
    res <- readGmt(file, uniqGenes=uniqGenes, namespace=namespace)
    return(res)
}
