\name{absLog10p}
\alias{absLog10p}
\title{
  Absolute base-10 logarithm of p-values
}
\description{
  The function returns the absolute values of base-10 logarithm of p-values.
}
\usage{
absLog10p(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector or matrix}
}
\details{
  The logarithm transformation of p-values is commonly used to visualize
  results from statistical tests. Although it may cause misunderstanding
  and therefore its use is disapproved by some experts, it helps to
  visualize and interpret results of statistical tests intuitively.

  The function transforms p-values with base-10 logarithm, and returns its
  absolute value. The choice of base 10 is driven by the simplicity of
  interpreting the results.
}
\value{
  Numeric vector or matrix.
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
testp <- runif(1000, 0, 1)
testp.al <- absLog10p(testp)

print(head(testp))
print(head(testp.al))
}
