\name{distance-methods}
\docType{methods}
\alias{distance}
\alias{distance-methods}
\alias{distance,Dna-method}
\title{Calculates absolute pairwise character difference matrix using a \code{Dna} object}
\description{
 Computes and returns an absolute pairwise character difference matrix from DNA sequences.
}

\usage{
\S4method{distance}{Dna}(x,subset=NULL,indels="sic")
}
\arguments{
 \item{x}{an object of class  \code{\link{Dna}}. }
 \item{subset}{a vector of integers in the range [1,nrow(x)], specifying which sequence(s) are used in the distance calculation. Only distance between selected sequence(s) and the rest of the sequences are calculated. If it is NULL, all comparisons are done. }
 \item{indels}{the indel coding method to be used. This must be one of "sic", "5th" or "missing". Any unambiguous substring can be given.  See also \sQuote{Details} }
}


\details{
  Available indel coding methods:
  \describe{
  	\item{\code{sic}:}{Treating gaps as a missing character and coding them separately following the simple indel coding method.}
    \item{\code{5th}:}{Treating gaps as a fifth state character.}
    \item{\code{missing}:}{Treating gaps as a missing character.}
}
}

\value{
  returns an object of class \code{dist}.
}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
Computes and returns an absolute pairwise character difference matrix from \code{Dna} objects.
}
}}


\references{
Giribet, G. and Wheeler, W.C. (1999) On gaps. \emph{Molecular Phylogenetics and Evolution} \bold{13}, 132-143.
  
Simmons, M., Ochoterena, H. (2000) Gaps as characters in sequence-based phylogenetic analyses. \emph{Systematic Biology} \bold{49}, 369-381.
		
  }

\seealso{
  \code{\link{indelcoder}} and \code{\link{subs}}
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\examples{
data("dna.obj")
x<-dna.obj[4:7,13:22,as.matrix=FALSE]

## Simple indel coding.
distance(x,indels="s")

## Gaps as 5th state characters.
distance(x,indels="5")

## Gaps as missing characters.
distance(x,indels="m")

\dontrun{
## Using 'subset'.
x<-dna.obj[4:10,13:22,as.matrix=FALSE]
distance(x, NULL)
distance(x, subset=c(1))
distance(x, subset=c(2,4))
}
}

\keyword{DNA ANALYSIS}



