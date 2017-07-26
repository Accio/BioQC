#' Convert a list of gene symbols into a gmtlist
#'
#' @param list: A named list with character vectors of genes. Names will become names of gene sets; character vectors will become genes
#' @param description: description; will be expanded to the same length of the list
#'
#' @examples
#' testVec <- list(GeneSet1=c("AKT1", "AKT2"),
#'                GeneSet2=c("MAPK1", "MAPK3"),
#'                GeneSet3=NULL)
#' testVecGmtlist <- as.gmtlist(testVec)
as.gmtlist <- function(list, description=NULL) {
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
                      list(name=names[i],
                           desc=desc,
                           genes=list[[i]])
                  })
    names(res) <- names
    return(GmtList(res))
}

readGmt <- function(filename) {
  stopifnot(file.exists(filename))
  lines <- readLines(filename)
  splitLines <- strsplit(lines, "\t")
  isValid <- sapply(splitLines, function(x) length(x)>=3)
  validLines <- splitLines[isValid]
  res <- lapply(validLines, function(x)  {
                    list(name=x[1], desc=x[2], genes=x[3:length(x)])
                })
  names(res) <- sapply(res, function(x) x$name)
  ## res <- .Call("read_gmt", filename, PACKAGE="BioQC")
  return(GmtList(res))
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
gmtlist2signedGenesets <- function(gmtlist, posPattern="_UP$", negPattern="_DN$",
                                   nomatch=c("ignore", "pos", "neg")) {
    nomatch <- match.arg(nomatch)
    genes <- lapply(gmtlist, function(x) x$genes)
    names <- sapply(gmtlist, function(x) x$name)
    isPos <- grepl(posPattern, names)
    isNeg <- grepl(negPattern, names)
    stemNames <- gsub(negPattern, "", gsub(posPattern, "", names))
    stemFactor <- factor(stemNames, levels=unique(stemNames))
    res <- tapply(seq(along=gmtlist), stemFactor, function(i) {
                      pos <- NULL
                      neg <- NULL
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
                      return(list(name=as.character(stemFactor)[i][1], pos=pos, neg=neg))
                  })
    names(res) <- levels(stemFactor)
    return(SignedGenesets(res))
}

#' Read signed GMT files
#'
#' @param filename A gmt file
#' @param posPattern: Pattern of positive gene sets
#' @param negPattern: Pattern of negative gene sets
#' @param nomatch: options to deal with gene sets that match to neither posPattern nor negPattern patterns
#'
#' @seealso \code{\link{gmtlist2signedGenesets}} for parameters \code{posPattern}, \code{negPattern}, and \code{nomatch}
#' @examples
#' testGmtFile <- system.file("extdata/test.gmt", package="BioQC")
#' testSignedGenesets.ignore <- readSignedGmt(testGmtFile, nomatch="ignore")
#' testSignedGenesets.pos <- readSignedGmt(testGmtFile, nomatch="pos")
#' testSignedGenesets.neg <- readSignedGmt(testGmtFile, nomatch="neg")
#' 
readSignedGmt <- function(filename, posPattern="_UP$", negPattern="_DN$",
                          nomatch=c("ignore", "pos", "neg")) {
    gmt <- readGmt(filename)
    res <- gmtlist2signedGenesets(gmt, posPattern=posPattern, negPattern=negPattern, nomatch=nomatch)
    return(res)
    
}
