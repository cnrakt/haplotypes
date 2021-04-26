\name{pairPhiST-methods}
\docType{methods}
\alias{pairPhiST}
\alias{pairPhiST-methods}
\alias{pairPhiST,dist-method}
\alias{pairPhiST,Dna-method}
\alias{pairPhiST,matrix-method}
\title{ Provides the pairwise PhiST between populations}
\description{
This function calculates pairwise PhiST between populations using the AMOVA framework.
}

\usage{
\S4method{pairPhiST}{Dna}(x,populations,indels="sic",nperm=99, negatives=FALSE, subset =NULL, 
showprogbar=TRUE)
\S4method{pairPhiST}{dist}(x,populations,nperm=99, negatives=FALSE, subset=NULL,showprogbar=TRUE)
\S4method{pairPhiST}{matrix}(x,populations,nperm=99,negatives=FALSE, subset=NULL,showprogbar=TRUE)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}, \code{"dist"} or \code{"matrix"}.}
 \item{populations}{a vector giving the populations, with one element per individual.}
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See \code{\link{distance}} for details. }
 \item{nperm}{the number of permutations. Set this to 0 to skip the permutation procedure.}
  \item{negatives}{boolean; if it is FALSE all negative PhiST values are replaced with zero. }
 \item{subset}{a vector of integers in the range [1, \code{length(unique(populations))}], only distances between selected population(s) and the rest of the populations are calculated. If it is NULL, all comparisons are done.}
 \item{showprogbar}{boolean; whether the progress bar is displayed or not  displayed.}

}

\details{
  The null distribution of pairwise PhiST under the hypothesis of no difference between the populations is obtained by permuting individuals between populations.   
}

\value{
	a list with following components: 
  	\item{\code{PhiST}}{a matrix giving the PhiST values between populations.}
	\item{\code{p}}{a matrix giving the p-values, or NULL if permutation test is not performed.}
}


\note{
	An internal code \code{Statphi} is taken from package  \pkg{ade4} version 1.7-8 without any modification, author Sandrine Pavoine. Function \code{amova} from package \pkg{pegas} is used internally to estimate variance components, author Emmanuel Paradis.

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
Excoffier, L., Smouse, P.E. and Quattro, J.M. (1992) Analysis of molecular variance inferred
from metric distances among DNA haplotypes: application to human mitochondrial DNA restriction
data. \emph{Genetics}, \bold{131}, 479-491.
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
pst<-pairPhiST(x, populations, nperm=0)
pst

##allow negative PhiST values
pst<-pairPhiST(x, populations, nperm=0, negatives=TRUE)
pst

##Gaps as missing characters.
pst<-pairPhiST(x, populations, indels="m", nperm=0, negatives=TRUE)
pst

##using subset, second population against others
pst <-pairPhiST(x, populations, nperm=0,subset=c(2))
pst


\dontrun{
## 999 permutations.
pst<-pairPhiST(x, populations, nperm=999,showprogbar=TRUE)
pst

## random populations
x<-dna.obj
populations<-sample(1:4,nrow(x),replace=TRUE)
pst<-pairPhiST(x, populations, nperm=999,showprogbar=TRUE)
pst

## populations based on clusters
x<-dna.obj
d<-distance(x)
hc<-hclust(d,method="ward.D")
populations<-cutree(hc,4)
pst<-pairPhiST(x, populations, nperm=999,showprogbar=TRUE)
pst
}

### Method for signature 'dist'.
x<-dna.obj
x<-dna.obj[c(1,20,21,26,27,28,30,3,4,7,13,14,15,16,23,24,25),,as.matrix=FALSE]
populations<-c("pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop2",
"pop2","pop2","pop2","pop2","pop3","pop4","pop3","pop4","pop4")
d<-distance(x)
pst<-pairPhiST(d, populations, nperm=0)
pst

### Method for signature 'matrix'.
x<-dna.obj
x<-dna.obj[c(1,20,21,26,27,28,30,3,4,7,13,14,15,16,23,24,25),,as.matrix=FALSE]
populations<-c("pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop2",
"pop2","pop2","pop2","pop2","pop3","pop4","pop3","pop4","pop4")
d<-as.matrix(distance(x))
pst<-pairPhiST(d, populations, nperm=0)
pst

}


\keyword{DNA}
\keyword{POPULATIONS}
