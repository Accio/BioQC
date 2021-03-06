% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readCurrentSignatures.R
\name{readCurrentSignatures}
\alias{readCurrentSignatures}
\title{Load current BioQC signatures}
\usage{
readCurrentSignatures(uniqGenes = TRUE, namespace = NULL)
}
\arguments{
\item{uniqGenes}{Logical, whether duplicated genes should be removed, passed
to \code{\link{readGmt}}.}

\item{namespace}{Character, namespace of the gene-set, or code{NULL}, passed
to \code{\link{readGmt}}}
}
\value{
A GmtList
}
\description{
Load current BioQC signatures
}
\examples{
readCurrentSignatures()
}
\seealso{
\code{\link{readGmt}}
}
