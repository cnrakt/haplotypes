\name{as.list-methods}
\docType{methods}

\alias{as.list-methods}
\alias{as.list,Dna-method}
\alias{as.list,Haplotype-method}
\alias{as.list,Parsimnet-method}
\title{Methods for function \code{as.list} in the Package \pkg{haplotypes}}
\description{
 Coerces an object to a list.
}

\usage{
\S4method{as.list}{Dna}(x)
\S4method{as.list}{Haplotype}(x)
\S4method{as.list}{Parsimnet}(x)
}

\arguments{
 \item{x}{ an object of class \code{\link{Dna}},\code{\link{Haplotype}} or \code{\link{Parsimnet}}.}
}

\details{
	If \code{x} is a \code{\link{Dna}} object, elements of the list are character vectors that contains the DNA sequences of length equal to corresponding value in the slot \code{seqlengths}. If \code{x} is \code{\link{Haplotype}} or \code{\link{Parsimnet}} objects, slots are converted to list elements.

		
	 }


\value{
   returns a list.
  }

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
coerces an object of class Dna to a list. 
}

\item{\code{signature(x = "Haplotype")}}{
coerces an object of class Haplotype to a list. 
}

\item{\code{signature(x = "Parsimnet")}}{
coerces an object of class Parsimnet to a list. 


}
}}

\examples{

data("dna.obj")

## Coercing a 'Dna' object to a list.
x<-dna.obj[1:3,as.matrix=FALSE]
as.list(x)

\dontrun{
## Coercing a 'Haplotype' object to a list.
x<-dna.obj
h<-haplotype(x)
as.list(h)

## Coercing a 'Parsimnet' object to a list.
x<-dna.obj
p<-parsimnet(x)
as.list(p)
}


}



\keyword{DNA ANALYSIS}
\keyword{HAPLOTYPE ANALYSIS}
\keyword{STATISTICAL PARSIMONY}

