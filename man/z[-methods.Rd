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
\S4method{[}{Dna}(x, i, j,..., drop = FALSE)
\S4method{[}{Dna}(x, i, j)<- value
}

\arguments{
 \item{x}{an object of class  \code{\link{Dna}} }
 \item{i, j}{elements to extract or replace.}
 \item{...}{
Additional arguments. In this context, ... is used primarily for as.matrix, which is a boolean. If as.matrix is TRUE (default), the function returns a matrix. If drop is also TRUE, and the subset has either a single row or column, the function will return a vector instead. If as.matrix is FALSE, the function returns an object of class \code{\link{Dna}}. The as.matrix argument should be specified as as.matrix=TRUE or as.matrix=FALSE within the call. If it is not specified, the function defaults to TRUE.
}
 \item{drop}{boolean; if TRUE and a single row or column is selected, the function returns a vector instead of a matrix. This is only applicable when as.matrix is TRUE.}
 \item{value}{a character vector or a character matrix. }
}

\details{
 The S4 method dispatch mechanism matches arguments based on those specified in the signature of the corresponding generic function. However, for some generics that include '...' in their signature, additional arguments can be incorporated in specific methods. Notably, the '[' function does not follow this pattern and restricts the arguments to those defined in its signature.
 In this context, the 'as.matrix' argument is not in the signature of the generic '[', so it is included within '...'. Then, within the body of the function, we check whether 'as.matrix' has been provided in the actual arguments when the function is called. If 'as.matrix' is not specified, the function defaults to 'TRUE', preserving the behavior of previous versions of the method.   
 }



\value{
 returns an object of class matrix, vector or  \code{\link{Dna}}.
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

# a vector object
x[1,1:4,drop=TRUE]

## Replace parts.
#"G" "C"
x[1,1:2] 	
x[1,1:2]<-c("A","T")
x[1,1:2]

}

\keyword{DNA ANALYSIS}
