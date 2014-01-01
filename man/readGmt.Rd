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
Gene set list
}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="ribiosBioQC")
gmt_list <- readGmt(gmt_file)
}
