\name{tfidf}
\alias{tfidf}
\alias{gsTfIdf}
\title{
Calculate TF-IDF
}
\description{
TF-IDF is a measure of informativeness of genes contributing to the gene set. The function calculates this measure.
}
\usage{
tfidf(x)
gsTfIdf(x)
}
\arguments{
  \item{x}{Any data structure that can be coereced into a matrix. Rows are genes (words) and columns are gene sets (documents).}
}
\value{
  \code{tfidf} returns a numeric vector. \code{gsTfIdf} returns a \code{\link{gslist-class}} object filled with weights.
}

\author{
Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
test.matrix <- matrix(c(1,0,0,1,1,1,1,1,0,0,0,1), nrow=3)
tfidf(test.matrix)
}
