\name{basecomp-methods}
\docType{methods}
\alias{basecomp}
\alias{basecomp-methods}
\alias{basecomp,Dna-method}
\title{ Calculates base composition }
\description{
 	Calculates base composition of \code{\link{Dna}} object.
}

\usage{
\S4method{basecomp}{Dna}(x)

}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}. }
 
}

\value{
  a matrix with sequence as rows, DNA bases as columns and frequencies as entries.

  }


\section{Methods}{
\describe{


\item{\code{signature(x = "Dna")}}{
	calculates base composition of \code{\link{Dna}} object.
}
}}


\examples{
data("dna.obj")
x <-dna.obj

## Calculating base compositions.
basecomp(x)



}

\keyword{DNA ANALYSIS}
