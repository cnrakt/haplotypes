\name{subs-methods}
\docType{methods}
\alias{subs}
\alias{subs-methods}
\alias{subs,Dna-method}
\title{Displays base substitutions}
\description{
This function displays all base substitutions. If \code{fifth=="TRUE"}, each gap is treated as a fifth state character.
}

\usage{
\S4method{subs}{Dna}(x,fifth=FALSE)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}. }
 \item{fifth}{ boolean; should gaps be treated as a fifth state character?}

}

\value{
	a list with three components: 
	 \describe{
  	\item{\code{subsmat}:}{a sequence matrix showing substitutions.}
	\item{\code{subs}:}{a list of matrices of the substitutions.}
	\item{\code{subsmnum}:}{total number of substitutions.}
}}




\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
}
}}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}. 


}



\examples{
data("dna.obj")
x<-dna.obj

## Base substitutions.
subs(x)

## Gaps are treated as a fifth state character.
subs(x,fifth=TRUE)

}


\keyword{DNA ANALYSIS}

