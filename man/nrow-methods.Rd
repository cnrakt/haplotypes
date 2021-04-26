\name{nrow-methods}
\docType{methods}

\alias{nrow-methods}
\alias{nrow,Dna-method}
\title{Returns the number of DNA sequences}
\description{
\code{nrow} returns the number of rows present in a matrix.
}
\usage{
\S4method{nrow}{Dna}(x)
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
\code{nrow} returns the number of rows present in the sequence matrix (number of sequences).
}
}}


\examples{
data("dna.obj")
x <-dna.obj

## Giving the number of sequences.
nrow(x)

}

\keyword{DNA ANALYSIS}
