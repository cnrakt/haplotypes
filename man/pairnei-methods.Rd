\name{pairnei-methods}
\docType{methods}
\alias{pairnei}
\alias{pairnei-methods}
\alias{pairnei,dist-method}
\alias{pairnei,Dna-method}
\alias{pairnei,matrix-method}
\title{Provides the average number of pairwise Nei's (D) differences between populations}
\description{
 Function provides pairwise Nei's raw number of nucleotide differences between populations. 
 
}

\usage{
\S4method{pairnei}{Dna}(x,populations,indels="sic",nperm=99, subset=NULL,showprogbar=FALSE)
\S4method{pairnei}{dist}(x,populations,nperm=99, subset=NULL,showprogbar=FALSE)
\S4method{pairnei}{matrix}(x,populations,nperm=99, subset=NULL,showprogbar=FALSE)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}, \code{"dist"} or \code{"matrix"}.}
 \item{populations}{a vector giving the populations, with one element per individual.}
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See \code{\link{distance}} for details. }
 \item{nperm}{the number of permutations. Set this to 0 to skip the permutation procedure.}
 \item{subset}{a vector of integers in the range [1, \code{length(unique(populations))}], only distances between selected population(s) and the rest of the populations are calculated. If it is NULL, all comparisons are done.}
 \item{showprogbar}{boolean; whether the progress bar is displayed or not  displayed.}

}

\details{
  The null distribution of pairwise Nei's differences under the hypothesis of no difference between the populations is obtained by permuting individuals between populations.  
}


\value{
	a list with following components: 
  	\item{\code{neidist}}{a matrix giving the average number of pairwise Nei's (D) differences between populations (below diagonal elements) and  average number of pairwise differences within populations (diagonal elements).}
	\item{\code{p}}{a matrix giving the p-values, or NULL if permutation test is not performed.}

}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
}

\item{\code{signature(x = "dist")}}{

}

\item{\code{signature(x = "matrix")}}{
}
}}

\references{
Nei, M. and Li, W. H. (1979) Mathematical model for studying genetic variation in terms of restriction endonucleases. \emph{Proceedings of the National Academy of Sciences of the United States of America } \bold{76}, 5269-5273.		
		
  }

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}



\examples{
data("dna.obj")


### Method for signature 'Dna'.
x<-dna.obj
x<-dna.obj[c(1,20,21,26,27,28,30,3,4,7,13,14,15,16,23,24,25),,as.matrix=FALSE]
populations<-c("pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop2",
"pop2","pop2","pop2","pop2","pop3","pop4","pop3","pop4","pop4")

##skip permutation testing
pn<-pairnei(x, populations, nperm=0)
pn

#Between populations
as.dist(pn$neidist)

#Within populations
diag(pn$neidist)

##Gaps as missing characters.
pn <-pairnei(x, populations, indels="m", nperm=0)
pn

##using subset, third population against others
pn<-pairnei(x, populations, nperm=0,subset=c(3))
pn

\dontrun{
## 999 permutations.
pn<-pairnei(x, populations, nperm=999, showprogbar=TRUE)
pn

## random populations
x<-dna.obj
populations<-sample(1:4,nrow(x),replace=TRUE)
pn<-pairnei(x, populations, nperm=999, showprogbar=TRUE)
pn

## populations based on clusters
x<-dna.obj
d<-distance(x)
hc<-hclust(d,method="ward.D")
populations<-cutree(hc,4)
pn<-pairnei(x, populations, nperm=999, showprogbar=TRUE)
pn
}

### Method for signature 'dist'.
x<-dna.obj
x<-dna.obj[c(1,20,21,26,27,28,30,3,4,7,13,14,15,16,23,24,25),,as.matrix=FALSE]
populations<-c("pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop2",
"pop2","pop2","pop2","pop2","pop3","pop4","pop3","pop4","pop4")
d<-distance(x)
pn<-pairnei(d, populations,nperm=0)
pn

### Method for signature 'matrix'.
x<-dna.obj
x<-dna.obj[c(1,20,21,26,27,28,30,3,4,7,13,14,15,16,23,24,25),,as.matrix=FALSE]
populations<-c("pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop2",
"pop2","pop2","pop2","pop2","pop3","pop4","pop3","pop4","pop4")
d<-as.matrix(distance(x))
pn<-pairnei(d, populations,nperm=0)
pn

}




\keyword{DNA}
\keyword{POPULATIONS}
