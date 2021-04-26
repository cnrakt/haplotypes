\name{parsimnet-methods}
\docType{methods}
\alias{parsimnet}
\alias{parsimnet-methods}
\alias{parsimnet,Dna-method}
\alias{parsimnet,dist-method}
\alias{parsimnet,matrix-method}
\title{Estimates gene genealogies using statistical parsimony}

\description{
 Function for estimating gene genealogies from DNA sequences or user provided absolute pairwise character difference matrix using statistical parsimony.
}

\usage{
\S4method{parsimnet}{Dna}(x,indels="sic",prob=.95)
\S4method{parsimnet}{dist}(x,seqlength,prob=.95)
\S4method{parsimnet}{matrix}(x,seqlength,prob=.95)

}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}, dist, or matrix.}
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See \code{\link{distance}} for details.}
  \item{seqlength}{an integer of length one giving the sequence length information (number of characters). }
  \item{prob}{a numeric vector of length one in the range [0.01, 0.99] giving the probability of parsimony as defined in Templeton et al. (1992). In order to set maximum connection steps to Inf (to connect all the haplotypes in a single network), set the probability to NULL. }
 
  }
  
   \details{
 The network estimation methods implemented in  \code{parsimnet} function finds one of the most parsimonious network (or sub-networks if connection between haplotypes exceeds the parsimony limit).  This is an implemetation of the TCS method proposed in Templeton et al. (1992) and Clement et al. (2002).\code{parsimnet} function generates an unambiguous haplotype network without loops. If more than one best networks found (results in ambiguous connections), only a network with the lowest average all-pairs distance is returned. Loops may occur only if they are present in initial haplotype distance matrix.
   }

\value{
  S4 methods for signature 'Dna', 'matrix' or 'dist' returns an object of class \code{\link{Parsimnet}}.}


\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
estimating gene genealogies from DNA sequences.
}

\item{\code{signature(x = "dist")}}{
estimating gene genealogies from distance matrix (dist object).
}

\item{\code{signature(x = "matrix")}}{
estimating gene genealogies from distance matrix.
}


}}



  \note{ 
  Duplicate names in the final distance matrices in slot d are renamed without warning. An internal function \code{.TempletonProb} is taken from package pegas version 0.6 without any modification, authors Emmanuel Paradis, Klaus Schliep. 
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}. 

}

\references{
	
	Clement, M., Q. Snell, P. Walker, D. Posada, and K. A. Crandall (2002) TCS: Estimating Gene Genealogies  \emph{in First IEEE International Workshop on High Performance Computational Biology (HiCOMB) } 
  
  Templeton, A. R., Crandall, K. A. and Sing, C. F. (1992) A cladistic
  analysis of phenotypic associations with haplotypes inferred from
  restriction endonuclease mapping and DNA sequence data. III. Cladogram
  estimation. \emph{Genetics}, \bold{132}, 619-635.

  
}




\seealso{
  \code{\link{network}}, \code{\link{plot-methods}} and \code{\link{pieplot-methods}}
  }


\examples{
	
\dontrun{
data("dna.obj")
x<-dna.obj

### Method for signature 'Dna'.
## statistical parsimony with 95% connection limit.
p<-parsimnet(x) 
p
plot(p)

## statistical parsimony with 99% connection limit. 
p<-parsimnet(x,prob=.99) 
p
# plot the first network
plot(p,net=1) 

## statistical parsimony with 99% connection limit.
#indels are coded as missing
p<-parsimnet(x,indels="m",prob=.99)
p
plot(p)


# statistical parsimony without connection limit.
p<-parsimnet(x,prob=NULL) 
p
plot(p)


# plot the first network
plot(p,net=1)

### Method for signature 'dist'.
d<-distance(x)
seqlength<-length(x)
## statistical parsimony with 95% connection limit for 113 character
p<-parsimnet(d,seqlength) 
p
plot(p)

### Method for signature 'matrix'.
d<-as.matrix(distance(x))
seqlength<-length(x)
## statistical parsimony with 95% connection limit for 113 character
p<-parsimnet(d,seqlength)
p
plot(p)

}
}


\keyword{STATISTICAL PARSIMONY}

