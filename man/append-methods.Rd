\name{append-methods}
\docType{methods}

\alias{append-methods}
\alias{append,Dna-method}
\title{Combines two \code{Dna} objects}
\description{
Combines two \code{\link{Dna}} objects.
}

\usage{
\S4method{append}{Dna}(x,values)

}
\arguments{
 \item{x}{an object of class \code{\link{Dna}}. }
 \item{values}{an object of class \code{\link{Dna}}. }

}

\value{
  an object of class \code{\link{Dna}}.
  }


\section{Methods}{
\describe{


\item{\code{signature(x = "Dna", values= "Dna")}}{
	combines two \code{Dna} objects.
}
}}


\examples{
	
data("dna.obj")

x<-dna.obj
y<-dna.obj
nrow(x)

## Combining two 'Dna' objects.
z<- append(x,y)
nrow(z)


}

\keyword{DNA ANALYSIS}
