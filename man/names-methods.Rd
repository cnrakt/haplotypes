\name{names-methods}
\docType{methods}
\alias{names-methods}
\alias{names<--methods}
\alias{names,Dna-method}
\alias{names<-,Dna-method}
\alias{names,Parsimnet-method}
\alias{names<-,Parsimnet-method}

\title{Function to get or set names of a \code{Dna} object or \code{Parsimnet} object}
\description{
 Function to get or set sequence names of \code{Dna} object or names of network in \code{Parsimnet} object.
}

\usage{
\S4method{names}{Dna}(x)
\S4method{names}{Parsimnet}(x)
\S4method{names}{Dna}(x)<-value
\S4method{names}{Parsimnet}(x)<-value
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}} or \code{\link{Parsimnet}}. }
 \item{value}{a character vector of the same length as number of sequence or networks.}
}

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
	 Function to get or set names of an object of Dna.
}
\item{\code{signature(x = "Parsimnet")}}{
	 Function to get or set names of networks in Parsimnet object.
}
}}


\examples{

data("dna.obj")

x<-dna.obj
x<-as.dna(x[1:4,1:6])

## Getting sequence names.
names(x)

## Setting sequence names.
names(x)<-c("u","v","z","y")
names(x)



x<-dna.obj
p<-parsimnet(x,prob=.99)

##Getting network names in parsimnet object
names(p)
## Setting network names names.
names(p)<-c("a","b","c","d","f","g")
names(p)

}


\keyword{DNA ANALYSIS}