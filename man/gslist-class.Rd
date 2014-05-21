\name{gslist-class}
\Rdversion{1.1}
\docType{class}
\alias{gslist-class}
\alias{coerce,data.frame,gslist-method}
\alias{coerce,gslist,data.frame-method}
\alias{coerce,gslist,matrix-method}
\alias{coerce,list,gslist-method}
\alias{gsClearWeights,gslist-method}
\alias{gsDesc,gslist-method}
\alias{gsGenes,gslist-method}
\alias{gsHasWeights,gslist-method}
\alias{gsIndices,gslist-method}
\alias{gsName,gslist-method}
\alias{gsUniqGenes,gslist-method}
\alias{gsWeights,gslist-method}
\alias{gsWeights<-,gslist,list-method}
\alias{show,gslist-method}

\title{Class \code{"gslist"}}
\description{
Data structure holding a collection (list) of \code{\link{gsitem-class}} objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gslist", ...)}. It is normally called internally by BioQC.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{List}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "data.frame", to = "gslist")} }
    \item{coerce}{\code{signature(from = "gslist", to = "data.frame")}}
    \item{coerce}{\code{signature(from = "gslist", to = "matrix")}}
    \item{coerce}{\code{signature(from = "list", to = "gslist")}}
    \item{gsClearWeights}{\code{signature(obj = "gslist")}}
    \item{gsDesc}{\code{signature(obj = "gslist")}}
    \item{gsGenes}{\code{signature(obj = "gslist")}}
    \item{gsHasWeights}{\code{signature(obj = "gslist")}}
    \item{gsIndices}{\code{signature(obj = "gslist")}}
    \item{gsName}{\code{signature(obj = "gslist")}}
    \item{gsUniqGenes}{\code{signature(obj = "gslist")}}
    \item{gsWeights}{\code{signature(obj = "gslist")}}
    \item{gsWeights<-}{\code{signature(obj = "gslist", value = "list")}}
    \item{show}{\code{signature(object = "gslist")}}
	 }
}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
showClass("gslist")
}
