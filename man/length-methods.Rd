\name{length-methods}
\docType{methods}

\alias{length-methods}
\alias{length,Dna-method}
\alias{length,Haplotype-method}
\alias{length,Parsimnet-method}
\title{Methods for function \code{length} in the package \pkg{haplotypes}}
\description{
Methods for function \code{length}.
}

\usage{
\S4method{length}{Dna}(x)
\S4method{length}{Haplotype}(x)
\S4method{length}{Parsimnet}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}, \code{\link{Haplotype}} or \code{\link{Parsimnet}}.}
}

\value{
   returns a non-negative integer vector.
  }

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
	returns the longest sequence length.
}

\item{\code{signature(x = "Haplotype")}}{
	 returns the number of haplotypes.
}

\item{\code{signature(x = "Parsimnet")}}{
	 returns the length of network(s). 
}
}}

\seealso{
  \code{\link{ncol-methods}}}


\examples{
data("dna.obj")
x<-dna.obj

## Longest sequence length
length(x)


## Total number of haplotypes
h<-haplotype(x)
length(h)

## Length of network(s)
p<-parsimnet(x,prob=.95)
# length of the network
length(p)
 
p<-parsimnet(x,prob=.99)
# length of the networks
length(p)

}



\keyword{DNA ANALYSIS}
\keyword{HAPLOTYPE ANALYSIS}
\keyword{STATISTICAL PARSIMONY}
