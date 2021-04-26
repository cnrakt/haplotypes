\name{boot.dna-methods}
\docType{methods}
\alias{boot.dna}
\alias{boot.dna-methods}
\alias{boot.dna,Dna-method}
\title{ Generates single bootstrap replicate}
\description{
 Methods for generating a single bootstrap replicate. }
 
 \usage{
\S4method{boot.dna}{Dna}(x,replacement=TRUE)
}

 \arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
 \item{replacement}{boolean; whether the sampling is done with replacement or without replacement.}
}
 
 \value{
  an object of class \code{\link{Dna}.}
  }


\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
generates single bootstrap replicate from a \code{Dna} object.
}
}}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}


\examples{
	
	
data("dna.obj")
x<-dna.obj

## Generating a bootstrap replicate.
# with replacement
bxr<-boot.dna(x)
image(bxr)

# without replacement
bx<-boot.dna(x,replacement=FALSE)
image(bx)
}

\keyword{DNA ANALYSIS}


