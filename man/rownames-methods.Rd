\name{rownames-methods}
\docType{methods}
\alias{rownames-methods}
\alias{rownames<--methods}
\alias{rownames,Dna-method}
\alias{rownames<-,Dna-method}
\alias{rownames,Parsimnet-method}
\alias{rownames<-,Parsimnet-method}

\title{Retrieve or set the row names}

\description{
 Function to get or set row names of a sequence matrix in a \code{Dna} object or distance matrix (or matrices) in a \code{Parsimnet} object.
}


\usage{
\S4method{rownames}{Dna}(x)
\S4method{rownames}{Parsimnet}(x)
\S4method{rownames}{Dna}(x)<-value
\S4method{rownames}{Parsimnet}(x)<-value
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}} or \code{\link{Parsimnet}}. }
 \item{value}{a character vector of the same length as number of sequences or a list of the same length as number of networks including vertex names for each network. See \sQuote{Examples} }
}



\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{

}

\item{\code{signature(x = "Parsimnet")}}{
}
}}



\examples{

data("dna.obj")
x<-dna.obj

### Method for signature 'Dna'.

## Getting sequence names.
rownames(x)

## Setting sequence names.
rownames(x)<-c(1:nrow(x))
rownames(x)

### Method for signature 'Parsimnet'.
x<-dna.obj

##single network
p<-parsimnet(x)

##Getting vertex names
rownames(p)

## Setting vertex names.
rownames(p)<-list(c(1:nrow(p@d[[1]])))
rownames(p)
plot(p)

## Multiple networks with 99% connection limit.
p<-parsimnet(x,prob=.99) 

## Getting vertex names
rownames(p)

## Setting vertex names.
rownames(p)<-list(1:9, 10, 11,12:13,14,15:16)
rownames(p)



}

\keyword{methods}

