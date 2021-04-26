\name{as.matrix-methods}
\docType{methods}

\alias{as.matrix-methods}
\alias{as.matrix,Dna-method}
\title{Methods for function \code{as.matrix} in the Package \pkg{haplotypes}}

\description{
 Coerces an object to a matrix.
}

\usage{
\S4method{as.matrix}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}

\value{
  returns a character matrix.
  }

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
coerces an object of class Dna to a matrix.
}
}}


\examples{
data("dna.obj")

## Coercing a 'Dna' object to a matrix.
x<-dna.obj[1:4,1:6,as.matrix=FALSE]
x
as.matrix(x)

\dontrun{
# gives the same result
dna.obj[1:4,1:6,as.matrix=TRUE] }

}


\keyword{DNA ANALYSIS}
