\name{as.DNAbin-methods}
\docType{methods}
\alias{as.DNAbin}
\alias{as.DNAbin-methods}
\alias{as.DNAbin,Dna-method}

\title{ Coerces an object to a  \code{DNAbin} object}
\description{

This function coerces \code{Dna} object to \code{\link{DNAbin}} \code{{ape}} object .
}

\usage{
\S4method{as.DNAbin}{Dna}(x, endgaps=TRUE)
}

\arguments{
  \item{x}{ an object of class \code{\link{Dna}}.}
  \item{endgaps}{ boolean; gaps at the end of the sequences are included if this is TRUE. }
}

\value{
  an object of class \code{\link{DNAbin}.}
  }

\section{Methods}{
\describe{
\item{\code{signature(x = "Dna")}}{
coerces a \code{\link{Dna}} object to a \code{DNAbin} object.
}

}}

\examples{

## Coercing a Dna object to a DNAbin object.
data("dna.obj")

x<-dna.obj
dBin<-as.DNAbin(x)
dBin

#gaps at the end removed
dBin<-as.DNAbin(x, endgaps=FALSE)
dBin


}
\keyword{DNA ANALYSIS}