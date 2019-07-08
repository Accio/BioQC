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
#' testVecGmtlist <- as.gmtlist(testVec)
#' 
#' @export
as.gmtlist <- function(list, description=NULL,
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

#' Read in gene-sets from a GMT file
#' 
#' @param filename Character, GMT file name
#' @param uniqGenes Logical, whether duplicated genes should be removed
#' @param namespace Character, namespace of the gene-set. It can be used to specify namespace or sources of the gene-sets. By default it is set as \code{NULL}, so no namespace is used and all gene-sets are assumed to come from the same unspecified namespace. The option can be helpful when gene-sets from multiple namespaces are jointly used.
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
#' @examples 
#' gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
#' gmt_list <- readGmt(gmt_file)
#' gmt_nonUniqGenes_list <- readGmt(gmt_file, uniqGenes=FALSE)
#' gmt_namespace_list <- readGmt(gmt_file, uniqGenes=FALSE, namespace="myNamespace")
#' 
#' @export
readGmt <- function(filename, uniqGenes=TRUE, namespace=NULL) {
  stopifnot(file.exists(filename))
  lines <- readLines(filename)
  splitLines <- strsplit(lines, "\t")
  isValid <- sapply(splitLines, function(x) length(x)>=3)
  validLines <- splitLines[isValid]
  res <- lapply(validLines, function(x)  {
    genes <- x[3:length(x)]
    if(uniqGenes) {
      genes <- unique(genes)
    }
    list(name=x[1], desc=x[2], genes=genes, namespace=namespace)
  })
  names(res) <- sapply(res, function(x) x$name)
  return(GmtList(res))
}


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


#' Convert gmtlist into a list of signed genesets
#'
#' @param gmtlist A gmtlist object, probably read-in by \code{readGmt}
#' @param posPattern Regular expression pattern of positive gene sets. It is trimmed from the original name to get the stem name of the gene set. See examples below.
#' @param negPattern Regular expression pattern of negative gene sets. It is trimmed from the original name to get the stem name of the gene set. See examples below.
#' @param nomatch Options to deal with gene sets that match neither positive nor negative patterns. ignore: they will be ignored (but not discarded, see details below); pos: they will be counted as positive signs; neg: they will be counted as negative signs
#' @return An S4-object of \code{SignedGenesets}, which is a list of signed_geneset, each being a two-item list; the first item is 'pos', containing a character vector of positive genes; and the second item is 'neg', containing a character vector of negative genes.
#'
#' Gene set names are detected whether they are positive or negative. If neither positive nor negative, nomatch will determine how will they be interpreted. In case of \code{pos} (or \code{neg}), such genesets will be treated as positive (or negative) gene sets.In case nomatch is set to \code{ignore}, the gene set will appear in the returned values with both positive and negative sets set to \code{NULL}.
#' 
#' @examples
#' testInputList <- list(list(name="GeneSetA_UP",genes=LETTERS[1:3]),
#'                list(name="GeneSetA_DN", genes=LETTERS[4:6]),
#'                list(name="GeneSetB", genes=LETTERS[2:4]),
#'                list(name="GeneSetC_DN", genes=LETTERS[1:3]),
#'                list(name="GeneSetD_UP", genes=LETTERS[1:3]))
#' testOutputList.ignore <- gmtlist2signedGenesets(testInputList, nomatch="ignore")
#' testOutputList.pos <- gmtlist2signedGenesets(testInputList, nomatch="pos")
#' testOutputList.neg <- gmtlist2signedGenesets(testInputList, nomatch="neg")
#' @export
gmtlist2signedGenesets <- function(gmtlist, posPattern="_UP$", negPattern="_DN$",
                                   nomatch=c("ignore", "pos", "neg")) {
    nomatch <- match.arg(nomatch)
    genes <- gsGenes(gmtlist)
    names <- gsName(gmtlist)
    hasNS <- hasNamespace(gmtlist)
    namespace <- gsNamespace(gmtlist)
    isPos <- grepl(posPattern, names)
    isNeg <- grepl(negPattern, names)
    stemNames <- gsub(negPattern, "", gsub(posPattern, "", names))
    stemFactor <- factor(stemNames, levels=unique(stemNames))
    res <- tapply(seq(along=gmtlist), stemFactor, function(i) {
                      pos <- NULL
                      neg <- NULL
                      if(hasNS) {
                        ns <- unique(namespace[i])
                        if(length(ns)>1) {
                          stop("Signed gene-sets found in multiple namespaces - please report the bug to the developer")
                        }
                      } else {
                        ns <- NULL 
                      }
                      if(!any(isPos[i]) && !any(isNeg[i])) {
                          otherGenes <- unique(unlist(genes[i]))
                          if (nomatch=="pos") {
                              pos <- otherGenes
                          } else if (nomatch=="neg") {
                              neg <- otherGenes
                          } else if (nomatch!="ignore") {
                              stop("should not be here")
                          }
                      } else {
                          if(any(isPos[i])) {
                              pos <- unique(unlist(genes[i][isPos[i]]))
                          }
                          if(any(isNeg[i])) {
                              neg <- unique(unlist(genes[i][isNeg[i]]))
                          }
                      }
                      return(list(name=as.character(stemFactor)[i][1], 
                                  pos=pos, 
                                  neg=neg,
                                  namespace=ns))
                  })
    names(res) <- levels(stemFactor)
    return(SignedGenesets(res))
}

#' Read signed GMT files
#'
#' @param filename A gmt file
#' @param posPattern Pattern of positive gene sets
#' @param negPattern Pattern of negative gene sets
#' @param nomatch options to deal with gene sets that match to neither posPattern nor negPattern patterns
#' @param uniqGenes Logical, whether genes should be made unique
#' @param namespace Character string or \code{NULL}, namespace of gene-sets
#' 
#' @seealso \code{\link{gmtlist2signedGenesets}} for parameters \code{posPattern}, \code{negPattern}, and \code{nomatch}
#' @examples
#' testGmtFile <- system.file("extdata/test.gmt", package="BioQC")
#' testSignedGenesets.ignore <- readSignedGmt(testGmtFile, nomatch="ignore")
#' testSignedGenesets.pos <- readSignedGmt(testGmtFile, nomatch="pos")
#' testSignedGenesets.neg <- readSignedGmt(testGmtFile, nomatch="neg")
#' @export
readSignedGmt <- function(filename, posPattern="_UP$", negPattern="_DN$",
                          nomatch=c("ignore", "pos", "neg"),
                          uniqGenes=TRUE, 
                          namespace=NULL) {
    gmt <- readGmt(filename, uniqGenes=uniqGenes, namespace=namespace)
    res <- gmtlist2signedGenesets(gmt, posPattern=posPattern, negPattern=negPattern, nomatch=nomatch)
    return(res)
    
}
