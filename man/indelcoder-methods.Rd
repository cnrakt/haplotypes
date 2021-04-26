\name{indelcoder-methods}
\docType{methods}
\alias{indelcoder}
\alias{indelcoder-methods}
\alias{indelcoder,Dna-method}
\title{Codes gaps}
\description{
Function for coding gaps separately. Only simple indel coding method is available in the current version.
}

\usage{
\S4method{indelcoder}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}


\value{
	a list with two components: 
	 \describe{
  	\item{\code{indels}:}{a matrix giving the indel positions (beginnings and ends) and lengths.}
	\item{\code{codematrix}:}{a binary matrix giving the indel codings. Missing values are denoted by -1.}
}}



\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
Function for coding gaps separately.

}
}}

\references{
Simmons, M., Ochoterena, H. (2000) Gaps as characters in sequence-based phylogenetic analyses. \emph{Systematic Biology} \bold{49}, 369-381.		
  }

\seealso{
  \code{\link{distance}}
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}


\examples{
data("dna.obj")
x<-dna.obj

## Simple indel coding.
indelcoder(x)

}

\keyword{DNA ANALYSIS}



