\name{as.networx-methods}
\docType{methods}
\alias{as.networx}
\alias{as.networx-methods}
\alias{as.networx,Parsimnet-method}
\title{ Coerces an object to a \code{networx} object}
\description{
This function coerces \code{Parsimnet} object to \code{\link{networx}} \code{{phangorn}} object .
}

\usage{
\S4method{as.networx}{Parsimnet}(x,net=1,...)
}

\arguments{
   \item{x}{an object of class \code{\link{Parsimnet}}.}
 	\item{net}{a numeric vector of length one indicating which network to convert.}
 	\item{...}{additional arguments to \code{\link{as.splits}}.}
}


\value{
  an object of class \code{\link{networx}.}
  }

\section{Methods}{
\describe{
\item{\code{signature(x = "Parsimnet")}}{
coerces a \code{\link{Parsimnet}} object to a \code{networx} object.
}

}}

\examples{

## Coercing a Parsimnet object to a networx object.
data("dna.obj")
x<-dna.obj
p<-parsimnet(x)
nx<-as.networx(p)
plot(nx, "2D")

}
\keyword{PHYLOGENETIC ANALYSIS}