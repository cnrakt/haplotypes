\name{homopoly-methods}
\docType{methods}
\alias{homopoly}
\alias{homopoly-methods}
\alias{homopoly,Dna-method}
\title{Provides the list of homoplastic indels and substitutions}
\description{
This function returs the list of homoplastic indels and substitutions.
}

\usage{
\S4method{homopoly}{Dna}(x,indels="sic",...)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
  \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See \code{\link{distance}} for details. }
\item{...}{additional arguments to \code{\link{parsimnet}}.}
  }


\value{
	a list with following components: 
  	\item{\code{indels}}{a character vector of homoplastic indels sitewise Consistency Index, names of the character vector gives the site of homoplastic indel.}
	\item{\code{subs}}{a character vector of homoplastic substitutions sitewise Consistency Index, names of the character vector gives the site of substitution.}
}




\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
}

}}



\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\examples{
data("dna.obj")

### Method for signature 'Dna'.
x<-dna.obj
homopoly(x)

}


\keyword{DNA}
\keyword{POPULATIONS}
