\name{tolower-methods, toupper-methods}
\docType{methods}
\alias{tolower-methods}
\alias{tolower,Dna-method}
\alias{toupper-methods}
\alias{toupper,Dna-method}
\title{ Convert sequence characters from upper to lower case or vice versa}
\description{
Convert sequence characters in a \code{Dna} object from upper to lower case or vice versa.
}

\usage{
\S4method{tolower}{Dna}(x)
\S4method{toupper}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}


\value{
	an object of class \code{\link{Dna}}.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{

}
}}
\examples{
	
## Coercing a list to a 'Dna' object.
seq1<-c("?","A","C","g","t","-","0","1")
seq2<-c("?","A","C","g","t","-","0","1","2")
seq3<-c("?","A","C","g","t","-","0","1","2","3")
x<-list(seq1=seq1,seq2=seq2,seq3=seq3) 
dna.obj<-as.dna(x)

#characters in Dna object
table(as.matrix(dna.obj))

##all lower case
lowc<-tolower(dna.obj)
#characters
table(as.matrix(lowc))

##all upper case
upc<-toupper(dna.obj)
#characters
table(as.matrix(upc))

}

\keyword{DNA ANALYSIS}
