\name{polymorp-methods}
\docType{methods}
\alias{polymorp}
\alias{polymorp-methods}
\alias{polymorp,Dna-method}
\title{Displays polymorphic sites (base substitutions and indels) between two sequences}
\description{
	This function displays the polymorphic sites (base substitutions and indels) between the two sequences.
}


\usage{
\S4method{polymorp}{Dna}(x,pair,indels="sic")
}
\arguments{
 \item{x}{an object of class  \code{\link{Dna}}. }
 \item{pair}{a vector of integers in the range [1,nrow(x)] of length two, specifying sequence pair. }
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See \code{\link{distance}} for details. }
}


\value{
	a list with two components: 
	 \describe{
  	\item{\code{indels}:}{a list of matrices of the indel regions if \code{indels=="sic"}. The component names of the list gives the position of the indels. }
	\item{\code{subst}:}{a list of matrices of the base substitutions. If \code{indels=="5th"}, each gap is treated as a base substitution. The component names of the list gives the position of the base substitutions. }
}}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
Showing base substitutions and indels between the two sequences.

}
}}

\seealso{
  \code{\link{indelcoder}} and \code{\link{subs}}
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\examples{
data("dna.obj")
x<-dna.obj

## Showing base substitutions and indels between seq1 and seq6.

# gaps are coded following the simple indel coding method
polymorp(x,c(1,6),indels="s")

# gaps are coded as a fifth state character
polymorp(x,c(1,6),indels="5")

# gaps are treated as missing character 
polymorp(x,c(1,6),indels="m")


}



\keyword{DNA ANALYSIS}
