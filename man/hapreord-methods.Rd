\name{hapreord-methods}
\docType{methods}
\alias{hapreord}
\alias{hapreord-methods}
\alias{hapreord,Haplotype-method}
\title{Reorders haplotypes according to the ordering factor}
\description{
 Reorders haplotypes according to the ordering factor.
}

\usage{
\S4method{hapreord}{Haplotype}(x,order=c(1:x@nhap))
}

\arguments{
 \item{x}{an object of class \code{\link{Haplotype}}. }
 \item{order}{a vector giving the order of haplotypes, with one element per haplotype. }
}

\value{
 returns an object of class \code{\link{Haplotype}}.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "Haplotype")}}{
	Reorders haplotypes.
}
}}

\seealso{
  \code{\link{haplotype}}}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\examples{
data("dna.obj")
x<-dna.obj[1:6,,as.matrix=FALSE]
# inferring haplotypes from DNA sequences
h<-haplotype(x)

## Reordering haplotypes.

# length of the argument 'order' is equal to the number of haplotypes
rh<-hapreord(h,order=c(4,3,1,2))
rh

}


\keyword{HAPLOTYPE ANALYSIS}