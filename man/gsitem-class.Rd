\name{gsitem-class}
\Rdversion{1.1}
\docType{class}
\alias{gsitem-class}
\alias{gsClearWeights,gsitem-method}
\alias{gsDesc,gsitem-method}
\alias{gsDesc<-,gsitem,character-method}
\alias{gsGenes,gsitem-method}
\alias{gsGenes<-,gsitem,character-method}
\alias{gsHasWeights,gsitem-method}
\alias{gsIndices,gsitem-method}
\alias{gsIndices<-,gsitem,integer-method}
\alias{gsName,gsitem-method}
\alias{gsName<-,gsitem,character-method}
\alias{gsWeights,gsitem-method}
\alias{gsWeights<-,gsitem,numeric-method}
\alias{length,gsitem-method}
\alias{show,gsitem-method}

\title{Class \code{"gsitem"}}
\description{
  Data structure to hold one gene set item
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gsitem", ...)}. Normally it is called by BioQC internally.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Character string of length 1, gene set name}
    \item{\code{desc}:}{Character string of length 1, gene set description}
    \item{\code{genes}:}{Vector of character strings, genes contained in the set}
    \item{\code{weights}:}{Vector of numeric values, relative weights of the genes. Its length must equal that of \code{genes}}
    \item{\code{indices}:}{Vector of integer values, indices of genes in a specific data set}
  }
}
\section{Methods}{
  \describe{
    \item{gsClearWeights}{\code{signature(obj = "gsitem")}: clear gene set weights}
    \item{gsDesc}{\code{signature(obj = "gsitem")}: get gene set description}
    \item{gsDesc<-}{\code{signature(obj = "gsitem", value = "character")}: set gene set description }
    \item{gsGenes}{\code{signature(obj = "gsitem")}: get genes}
    \item{gsGenes<-}{\code{signature(obj = "gsitem", value = "character")}: set genes }
    \item{gsHasWeights}{\code{signature(obj = "gsitem")}: logical, whether gene sets have weights }
    \item{gsIndices}{\code{signature(obj = "gsitem")}: get indices of genes }
    \item{gsIndices<-}{\code{signature(obj = "gsitem", value = "integer")}: set indices of genes }
    \item{gsName}{\code{signature(obj = "gsitem")}: get gene set name}
    \item{gsName<-}{\code{signature(obj = "gsitem", value = "character")}: set gene set name}
    \item{gsWeights}{\code{signature(obj = "gsitem")}: get gene weights }
    \item{gsWeights<-}{\code{signature(obj = "gsitem", value = "numeric")}: set gene weights }
    \item{length}{\code{signature(x = "gsitem")}: number of genes }
    \item{show}{\code{signature(object = "gsitem")}: print method }
	 }
}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com> 
}

\seealso{
\code{\link{gslist-class}}
}
\examples{
showClass("gsitem")
}
