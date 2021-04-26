\name{ncol-methods}
\docType{methods}

\alias{ncol-methods}
\alias{ncol,Dna-method}
\title{Returns the length of the longest DNA sequence}
\description{
\code{ncol} returns the number of columns present in a matrix.
}

\usage{
\S4method{ncol}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}

\value{
	an integer of length one.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
\code{ncol} returns the number of columns present in the sequence matrix (length of the longest DNA sequence).
}
}}

\seealso{
  \code{\link{length-methods}}}

\examples{
data("dna.obj")
x <-dna.obj

## Giving the length of the longest sequence.
ncol(x) 
# gives the same result
length(x) 
}

\keyword{DNA ANALYSIS}

