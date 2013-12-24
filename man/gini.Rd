\name{gini}
\alias{gini}
\title{
  Calculate Gini Index of a numeric vector
}
\description{
  Calculate the Gini index of a numeric veector
}
\usage{
gini(x, na.rm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector.}
  \item{na.rm}{Logical. If set to \code{TRUE}, \code{NA} values are omitted.}
}
\details{
  The Gini index (Gini coefficient) is a measure of statistical dispersion. A Gini
  coefficient of zero expresses perfect equality where all values are
  the same. A Gini coefficient of one expresses maximal inequality among values.
}
\value{
  A numeric value between 0 and 1.
}
\references{
  Gini. C. (1912) \emph{Variability and Mutability}, C. Cuppini, Bologna
  156 pages. 
}
\author{
  Jitao David Zhang
}
\examples{

testValues <- runif(100)
gini(testValues)
}
