\name{[-methods}
\docType{methods}
\alias{[-methods}
\alias{[,Dna-method}
\alias{[<-,Dna-method}
\alias{[,Dna,ANY,ANY,ANY-method}
\alias{[<-,Dna,ANY,ANY,ANY-method}
\title{Extract or replace parts of an object of class \code{Dna}}
\description{
Operators acting on sequence matrix to extract or
  replace parts.
  }

  \usage{
\S4method{[}{Dna}(x, i, j, as.matrix = TRUE)
\S4method{[}{Dna}(x, i, j)<- value
}

\arguments{
 \item{x}{an object of class  \code{\link{Dna}} }
 \item{i, j}{elements to extract or replace.}
 \item{as.matrix}{boolean; If TRUE, function returns a matrix, else, function returns an object of class  \code{\link{Dna}}.}
 \item{value}{a character vector or a character matrix. }
 
 
}

\value{
 returns an object of class matrix or  \code{\link{Dna}}.

}
  
\section{Methods}{
\describe{

\item{\code{signature(x = "Dna", i = "ANY", j = "ANY", drop = "ANY")}}{
}
}}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\seealso{
\code{\link{Dna}}
}
\examples{

data("dna.obj")
x<-dna.obj

## Extract parts.
# a matrix object
x[1:2,1:3]

# a Dna object, as.dna(x[1:2,1:3]) gives the same result
x[1:2,1:3,as.matrix=FALSE] 
  
## Replace parts.
#"G" "C"
x[1,1:2] 	
x[1,1:2]<-c("A","T")
x[1,1:2]



}


\keyword{DNA ANALYSIS}
