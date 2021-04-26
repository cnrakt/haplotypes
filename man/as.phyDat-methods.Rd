\name{as.phyDat-methods}
\docType{methods}
\alias{as.phyDat}
\alias{as.phyDat-methods}
\alias{as.phyDat,Dna-method}


\title{ Coerces an object to a \code{phyDat} object}
\description{

This function coerces \code{Dna} object to \code{\link{phyDat}} \code{{phangorn}} object .
}


\usage{
\S4method{as.phyDat}{Dna}(x, indels="sic",...)
}

\arguments{
 \item{x}{an object of class  \code{\link{Dna}}. }
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See also \sQuote{Details} }
   \item{...}{additional arguments to \code{\link{as.phyDat}}.}
  }
  
  \details{
  Available indel coding methods:
  \describe{
  	\item{\code{sic}:}{Treating gaps as a missing character and coding them separately following the simple indel coding method.}
    \item{\code{5th}:}{Treating gaps as a fifth state character.}
    \item{\code{missing}:}{Treating gaps as a missing character.}
}
}

  

\value{
  an object of class \code{\link{phyDat}.}
  }

\section{Methods}{
\describe{
\item{\code{signature(x = "Dna")}}{
coerces a \code{\link{Dna}} object to a \code{phyDat} object.
}

}}

\examples{


data("dna.obj")
x<-dna.obj


## Coercing a Dna object to a phyDat object.
# Simple indel coding.
phyd<-as.phyDat(x)
phyd

# Gaps as 5th state characters.
phyd<-as.phyDat(x,indels="5")
phyd

# Gaps as 5th state characters.
phyd<-as.phyDat(x,indels="m")
phyd


}
\keyword{DNA ANALYSIS}