\name{unique-methods}
\docType{methods}
\alias{unique-methods}
\alias{unique,Dna-method}
\alias{unique,Dna,ANY-method}

\title{Extract Unique Sequences}
\description{
 \code{unique} returns a list with duplicate sequences removed. 
}

\usage{
\S4method{unique}{Dna}(x,gaps=FALSE)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
 \item{gaps}{boolean; gaps are removed if this is FALSE.}

}

\details{
	
	This function behaves somehow similar to \code{\link{haplotype}}, however indels and missing characters are not taken into account.
  
}


\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{

}
}}


\examples{
	
data("dna.obj")
x<-dna.obj[1:6,,as.matrix=FALSE]

##gaps removed.
unique(x)

##gaps not removed.
unique(x,gaps=TRUE)

##unique vs. haplotype
#unique returns 5 unique sequences.
unique(x)
length(unique(x))

#haplotype returns 4 unique haplotypes with simple indel coding.
h<-haplotype(x)
as.list(as.dna(h))
length(h)

#haplotype returns 3 unique haplotypes with gaps as missing.
h<-haplotype(x, indels="m")
as.list(as.dna(h))
length(h)

}

\keyword{DNA ANALYSIS}
