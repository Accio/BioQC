\name{readGmt}
\alias{readGmt}
\title{Read in gene sets from a GMT file}
\description{
Read in gene sets from a GMT file
}
\usage{
readGmt(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{GMT file name}
}
\value{
A gene set list, wrapped in a S3-class \code{gmtlist}. Each list item is
a list with three items: gene set name (\code{name}), gene set
description (\code{desc}), and gene list (a character vector,
\code{genes}). 
}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt_list <- readGmt(gmt_file)
}
