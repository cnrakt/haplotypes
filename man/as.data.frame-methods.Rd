\name{as.data.frame-methods}
\docType{methods}

\alias{as.data.frame-methods}
\alias{as.data.frame,Dna-method}
\title{Coerces a \code{Dna} object to a data.frame}
\description{
Coerces an object to a data.frame.
}
\usage{
\S4method{as.data.frame}{Dna}(x)
}

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
}

\value{
   returns a data frame.
  }


\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
coerces a Dna object to a data.frame. }
	}
}



\examples{

data("dna.obj")

x<-dna.obj
x<-as.dna(x[1:4,1:6])

## Coercing a 'Dna' object to a data.frame.
df<-as.data.frame(x)
df

# TRUE
is.data.frame(df) 

\dontrun{
# gives the same result 
df<-as.data.frame(x@sequence) 
df
}
}


\keyword{DNA ANALYSIS}




