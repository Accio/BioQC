\name{entropy}
\alias{entropy}
\alias{entropyDiversity}
\alias{entropySpecificity}
\alias{sampleSpecialization}
\title{
  Shannon entropy and related concepts
}
\description{
  These functions calculate Shannon entropy and related concepts,
  including diversity, specificity, and specialization. They can be used
  to quantify gene expression profiles.
}
\usage{
entropy(vector)
entropyDiversity(mat, norm=FALSE)
entropySpecificity(mat, norm=FALSE)
sampleSpecialization(mat, norm=TRUE)
}
\arguments{
  \item{vector}{A vector of numbers, or characters. Discrete probability
    of each item is calculated and the Shannon entropy is returned.}
  \item{mat}{A matrix (usually an expression matrix), with genes
    (features) in rows and samples in columns.}
  \item{norm}{Logical value. If set to \code{TRUE} the scores will be
    normalized between 0 and 1.}
}
\details{
  Shannon entropy can be used as measures of gene expression
  specificity, as well as measures of tissue diversity and
  specialization. See references below.

  We use \code{2} as base for the entropy calculation, because in this
  base the unit of entropy is \emph{bit}.
}
\value{
  \code{entropy} returns one entropy value. \code{entropyDiversity} and
  \code{sampleSpecialization} returns a vector as long as the column
  number of the input matrix. \code{entropySpecificity} returns a vector
  of the length of the row number of the input matrix, namely the
  specificity score of genes.
}
\references{
  Martinez and Reyes-Valdes (2008) Defining diversity, specialization, and
  gene specificity in transcriptomes through information theory. PNAS
  105(28):9709--9714
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
myVec0 <- 1:9
entropy(myVec0) ## log2(9)
myVec1 <- rep(1, 9)
entropy(myVec1)

myMat <- rbind(c(3,4,5),c(6,6,6), c(0,2,4))
entropySpecificity(myMat)
entropySpecificity(myMat, norm=TRUE)
entropyDiversity(myMat)
entropyDiversity(myMat, norm=TRUE)
sampleSpecialization(myMat)
sampleSpecialization(myMat,norm=TRUE)

myRandomMat <- matrix(runif(1000), ncol=20)
entropySpecificity(myRandomMat)
entropySpecificity(myRandomMat, norm=TRUE)
entropyDiversity(myRandomMat)
entropyDiversity(myRandomMat, norm=TRUE)
sampleSpecialization(myRandomMat)
sampleSpecialization(myRandomMat,norm=TRUE)

}
