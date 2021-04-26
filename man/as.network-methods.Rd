\name{as.network-methods}
\docType{methods}
\alias{as.network}
\alias{as.network-methods}
\alias{as.network,Parsimnet-method}
\title{ Coerces an object to a \code{network} object}
\description{
This function coerces \code{Parsimnet} object to \code{\link{network}} \code{{network}} object .
}


\usage{
\S4method{as.network}{Parsimnet}(x,net=1,...)
}

\arguments{
   \item{x}{an object of class \code{\link{Parsimnet}}.}
 	\item{net}{a numeric vector of length one indicating which network to convert.}
 	\item{...}{additional arguments to function \code{\link{network}}.}
}

\value{
  an object of class \code{\link{network}.}
  }

\section{Methods}{
\describe{
\item{\code{signature(x = "Parsimnet")}}{
coerces a \code{\link{Parsimnet}} object to a \code{network} object.
}

}}

\examples{

## Coercing a Parsimnet object to a network object.
data("dna.obj")
x<-dna.obj
p<-parsimnet(x)
n<-as.network(p)

#Fourth network (with only two edges)
p<-parsimnet(x,prob=.99)
n<-as.network(p,net=4)



}
\keyword{PHYLOGENETIC ANALYSIS}