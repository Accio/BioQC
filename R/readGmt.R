readSingleGmt <- function(filename, uniqGenes=TRUE, namespace=NULL) {
  if(!file.exists(filename))
    stop(sprintf("File %s does not exist", filename))
  lines <- readLines(filename)
  splitLines <- strsplit(lines, "\t")
  isValid <- sapply(splitLines, function(x) length(x)>=3)
  validLines <- splitLines[isValid]
  res <- lapply(validLines, function(x)  {
    genes <- x[3:length(x)]
    genes <- genes[!genes %in% ""]
    if(uniqGenes) {
      genes <- unique(genes)
    }
    list(name=x[1], desc=x[2], genes=genes, namespace=namespace)
  })
  names(res) <- sapply(res, function(x) x$name)
  return(res)
}

#' Read in gene-sets from a GMT file
#' 
#' @param ... Named or unnamed characater string vector, giving file names of one or more GMT format files. 
#' @param uniqGenes Logical, whether duplicated genes should be removed
#' @param namespace Character, namespace of the gene-set. It can be used to specify namespace or sources of the gene-sets. If  \code{NULL} is given, so no namespace is used and all gene-sets are assumed to come from the same unspecified namespace. The option can be helpful when gene-sets from multiple namespaces are jointly used.
#' 
#' @return A \code{GmtList} object, which is a S4-class wrapper of a list. Each 
#' element in the object is a list of (at least) three items: 
#' \itemize{
#'   \item gene-set name (field \code{name}), character string, accessible with \code{\link{gsName}}
#'   \item gene-set description (field \code{desc}), character string, accessible with \code{\link{gsDesc}}
#'   \item genes (field \code{genes}), a vector of character strings, , accessible with \code{\link{gsGenes}}
#'   \item namespace (field \code{namespace}), accessible with \code{\link{gsNamespace}}
#' }
#' 
#' @note 
#' Currently, when \code{namespace} is set as \code{NULL}, no namespace is used. This may change in the future, since we may use file base name as the default namespace.
#' 
#' @examples 
#' gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
#' gmt_list <- readGmt(gmt_file)
#' gmt_nonUniqGenes_list <- readGmt(gmt_file, uniqGenes=FALSE)
#' gmt_namespace_list <- readGmt(gmt_file, uniqGenes=FALSE, namespace="myNamespace")
#' 
#' ## suppose we have two lists of gene-sets to read in
#' test_gmt_file <- system.file("extdata/test.gmt", package="BioQC")
#' gmt_twons_list <- readGmt(gmt_file, test_gmt_file, namespace=c("BioQC", "test"))
#' ## alternatively
#' gmt_twons_list <- readGmt(BioQC=gmt_file, test=test_gmt_file)
#' 
#' @export
readGmt <- function(..., uniqGenes=TRUE, namespace=NULL) {
  files.list <- list(...)
  files <- unlist(files.list, use.names=TRUE)
  fnames <- names(files)
  if(!is.null(fnames)) {
    if(!is.null(namespace))
      warning("'file' have names - the namespace option is ignored")
    namespace <- fnames
  } else if(!is.null(namespace)) {
    namespace <- rep(namespace, length.out=length(files))
  }
  gsList <- lapply(seq(along=files), function(i) {
    readSingleGmt(files[i], uniqGenes=uniqGenes, namespace = namespace[i])
  })
  res <- GmtList(unlist(gsList, use.names=TRUE, recursive=FALSE))
  return(res)
}
