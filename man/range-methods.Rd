\name{range-methods}
\docType{methods}
\alias{range-methods}
\alias{range,Dna-method}
\title{Returns the minimum and maximum lengths of the DNA sequences}
\description{
\code{range} returns the lengths of shortest and longest DNA sequences.
}

\usage{
\S4method{range}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}

\value{
	an integer of length two.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
\code{range} 
}
}}

\seealso{
  \code{\link{length-methods}}}

\examples{
data("dna.obj")
x <-dna.obj

## shortest and longest DNA sequence lengths
range(x) 

}

\keyword{DNA ANALYSIS}

