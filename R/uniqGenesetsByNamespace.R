#' Make names of gene-sets unique by namespace, and member genes of gene-sets unique
#' 
#' @param gmtList A \code{GmtList} object, probably from \code{\link[BioQC]{readGmt}}. The object must have namespaces defined by \code{setNamespace}.
#' 
#' The function make sure that
#' \itemize{
#'   \item names of gene-sets within each namespace are unique, by merging gene-sets with duplicated names
#'   \item genes within each gene-set are unique, by removing duplicated genes
#' }
#' 
#' Gene-sets with duplicated names and different \code{desc} are merged, \code{desc} are made unique, and in case of multiple values, concatenated (with \code{|} as the collapse character).
#' 
#' @return A \code{GmtList} object, with unique gene-sets and unique gene lists. If not already present, a new item \code{namespace} is appended to each \code{list} element in the \code{GmtList} object, recording the namespace used to make gene-sets unique. The order of the returned \code{GmtList} object is given by the unique gene-set name of the input object.
#' 
#' @importFrom methods is
#' 
#' @examples 
#' myGmtList <- GmtList(list(list(name="GeneSet1", desc="Namespace1", genes=LETTERS[1:3]),
#'   list(name="GeneSet2", desc="Namespace1", genes=rep(LETTERS[4:6],2)),
#'   list(name="GeneSet1", desc="Namespace1", genes=LETTERS[4:6]),
#'   list(name="GeneSet3", desc="Namespace2", genes=LETTERS[1:5])))
#'  
#' print(myGmtList)
#' myGmtList <- setNamespace(myGmtList, namespace=function(x) x$desc)
#' myUniqGmtList <- uniqGenesetsByNamespace(myGmtList)
#' print(myUniqGmtList)
#' 
#' @export
uniqGenesetsByNamespace <- function(gmtList) {
  stopifnot(is(gmtList, "GmtList"))
  if(!hasNamespace(gmtList)) {
    stop(paste0("Gene-set namespace is not defined. ",
                "Please run 'setNamespace' or 'setDescAsNamespace' first to define namespace."))
  }
  
  inputNamespace <- gsNamespace(gmtList)
  cg <- factor(inputNamespace)
  
  inputNames <- gsName(gmtList)
  inputDescs <- gsDesc(gmtList)
  geneList <- gsGenes(gmtList)
  gsSize <- sapply(geneList, length)
  cgFac <- rep(cg, gsSize)
  tbl <- data.frame(namespace=cgFac,
                    name=rep(gsName(gmtList), gsSize),
                    genes=as.character(unlist(geneList)),
                    stringsAsFactors = FALSE)
  nameFactor <- factor(tbl$name)
  tbl$nameFac <- as.integer(nameFactor)
  tbl$cgFac <- as.integer(cgFac)
  tbl <- unique(tbl) ## unique combination of namespace, name, genes
  splitGenes <- with(tbl,
                     split(genes,
                           list(cgFac, nameFac), 
                           drop=TRUE))
  resList <- lapply(seq(along=splitGenes),
                    function(i) {
                      genes <- splitGenes[[i]]
                      intName <- names(splitGenes)[i]
                      intNameSplit <- strsplit(intName, "\\.")[[1]]
                      cgLevel <- as.integer(intNameSplit[1])
                      nameLevel <- as.integer(intNameSplit[2])
                      resNamespace <- levels(cg)[cgLevel]
                      resName <- levels(nameFactor)[nameLevel]
                      matchingDesc <- inputDescs[inputNames==resName &
                                                   inputNamespace==resNamespace]
                      resDesc <- paste(unique(matchingDesc),
                                       collapse="|")
                      res <- list(name=resName, 
                                  desc=resDesc, 
                                  genes=genes,
                                  namespace=resNamespace)
                      return(res)
                    })
  names(resList) <- sapply(resList, function(x) x$name)
  oResList <- resList[unique(inputNames)]
  res <- GmtList(oResList)
  return(res)
}