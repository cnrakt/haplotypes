\name{as.numeric-methods}
\docType{methods}

\alias{as.numeric-methods}
\alias{as.numeric,Dna-method}
\title{Coerces a \code{Dna} object to a numeric matrix}
\description{Converts a character matrix to a numeric matrix.
}

\usage{
\S4method{as.numeric}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}


\details{
 Function \code{as.numeric()} coerces the character matrix in the slot \code{sequence} to a numeric matrix. 
Lower or upper case characters "?","A","C","G","T","-" are converted to integers 0,1,2,3,4,5, respectively. 
}

\value{
 returns a numeric matrix.
  }

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
coerces a Dna object to a numeric matrix.
}
}}


\examples{
x<-matrix(c("?","A","C","g","t","-","0","1","2","3","4","5"),4,6)
rownames(x)<-c("seq1","seq2","seq3","seq4") 
x<-as.dna(x)
# original character matrix 
as.matrix(x) 

## Coercing a 'Dna' object to a numeric matrix.
# numeric matrix
as.numeric(x)

}

\keyword{DNA ANALYSIS}



