\name{haplotype-methods}
\docType{methods}
\alias{haplotype}
\alias{haplotype-methods}
\alias{haplotype,Dna-method}
\alias{haplotype,dist-method}
\alias{haplotype,matrix-method}
\title{Methods for function \code{haplotype} in the package \pkg{haplotypes}}
\description{
Collapses identical DNA sequences into haplotypes or infering haplotypes using user provided absolute pairwise character difference matrix. 
}

\usage{
\S4method{haplotype}{Dna}(x,indels="sic")
\S4method{haplotype}{dist}(x)
\S4method{haplotype}{matrix}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}, \code{dist}, or \code{matrix}.}
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See \code{\link{distance}} for details. }
}


\value{
  \code{haplotype} returns an object of class \code{\link{Haplotype}}, \code{\link{as.list-methods}} can be used to coerce the object to a list.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
Inferring haplotypes from DNA sequences.
}

\item{\code{signature(x = "dist")}}{
Inferring haplotypes using an absolute pairwise character difference matrix (dist object).
}

\item{\code{signature(x = "matrix")}}{
Inferring haplotypes using an absolute pairwise character difference matrix.
}
}}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}


\examples{
	
data("dna.obj")
x<-dna.obj[1:6,,as.matrix=FALSE]

##Inferring haplotypes using 'Dna' object.
# coding gaps using simple indel coding method
h<-haplotype(x,indels="sic")
h

# giving DNA sequences of haplotypes
as.dna(h) 

\dontrun{
	
## Slots of an object Haplotype	
h@haplist  #haplotype list (names)
h@hapind  #haplotype list (index)
h@uniquehapind  #getting index of the first occurrence of haplotypes
h@sequence  #DNA sequences of haplotypes
h@d  #distance matrix of haplotypes
h@freq  #haplotype frequencies
h@nhap  #total number of haplotypes
}

## Inferring haplotypes using dist object. 
d<-distance(x)
h<-haplotype(d) 
h
\dontrun{
# returns an error message 
as.dna(h)
} 

## Inferring haplotypes using distance matrix.
d<-as.matrix(distance(x))
h<-haplotype(d) 
h
\dontrun{
# returns an error message
as.dna(h)
}  

}

\keyword{HAPLOTYPE ANALYSIS}