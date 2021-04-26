\name{remove.gaps-methods}
\docType{methods}
\alias{remove.gaps}
\alias{remove.gaps-methods}
\alias{remove.gaps,Dna-method}
\title{ Removing gaps from \code{Dna} object}
\description{
 Removing gaps("-") from \code{Dna} object
}
\usage{
\S4method{remove.gaps}{Dna}(x,entire.col=FALSE)
}


\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
 \item{entire.col}{boolean; entire columns with gaps are removed if this is TRUE. See also \sQuote{Details}. }

}

\details{
	
	If entire.col==TRUE, alignment is preserved. If it is FALSE, end gaps are introduced to sequence matrix. 
  
}


\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{

}
}}

\value{
	an object of class  \code{\link{Dna}}.

}


\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\examples{

data("dna.obj")

## original data
x<-dna.obj 
range(x)
x@seqlengths

## Only gaps '-' are removed from sequences.
x<-remove.gaps(dna.obj, entire.col=FALSE)
range(x)
x@seqlengths

## entire columns with gaps are removed.
x<-remove.gaps(dna.obj, entire.col=TRUE)
range(x)
x@seqlengths


}
