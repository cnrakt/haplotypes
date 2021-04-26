\name{grouping-methods}
\docType{methods}
\alias{grouping}
\alias{grouping-methods}
\alias{grouping,Haplotype-method}
\title{Groups haplotypes according to the grouping variable (populations, species, etc.)
	}
\description{
Function for creating a matrix with haplotypes as rows, grouping factor (populations, species, etc.) as columns and abundance as entries.
}

\usage{
\S4method{grouping}{Haplotype}(x,factors)
}

\arguments{
 \item{x}{an object of class  \code{\link{Haplotype}}. }
 \item{factors}{a vector or factor giving the grouping variable (populations, species, etc.), with one element per individual. }
}

\value{
	a list with two components: 
	 \describe{
  	\item{\code{hapmat}:}{a matrix with haplotypes as rows, levels of the grouping factor (populations, species, etc.) as columns and abundance as entries.}

	\item{\code{hapvec}:}{a vector giving the haplotype identities of individuals. }
}}


\section{Methods}{
\describe{

\item{\code{signature(x = "Haplotype")}}{
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

## Grouping haplotypes.
# character vector 'populations' is a grouping factor.
populations<-c("pop1","pop1","pop2","pop3","pop3","pop3") 

# length of the argument 'factor' is equal to the number of sequences
g<-grouping(h,factors=populations) 
g

}

\keyword{HAPLOTYPE ANALYSIS}

