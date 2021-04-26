\name{image-methods}
\docType{methods}
\alias{image-methods}
\alias{image,Dna-method}
\title{ Display DNA Sequence }
\description{
 Display an image of DNA sequences .

}



\usage{
\S4method{image}{Dna}(x,all=FALSE,fifth=TRUE,
col=c("#BFBFBF","#0B99FD","#FD0B0B","#11A808","#F5FD0B","#F8F8FF"),
chars=TRUE,cex=1,show.names=TRUE,show.sites=TRUE,xlab="",ylab="",...)
 }

\arguments{
 \item{x}{an object of class \code{\link{Dna}}.}
  \item{all}{boolean; should entire sequence be displayed or only the polymorphic sites? }
  \item{fifth}{boolean; if all==FALSE, should gaps be displayed? }
	\item{col}{ an integer or character vector for the colors. By default it is blue for "A", red for "C", green for "G", yellow for "T", white for "-",  and grey for "?".}
	\item{chars}{boolean; should characters be displayed on image?}
	\item{cex}{a numeric vector of expansion factor for characters. }	
 	\item{show.names}{boolean; should sequence names be displayed on the left side.}
 	\item{show.sites}{boolean; should site labels be displayed on the bottom side.}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{...}{additional arguments to \code{\link{image.default}}.}
  }

\section{Methods}{
\describe{

\item{\code{signature(x = "Dna")}}{
	Display an image of Dna objects
}
}}

\seealso{
 \code{\link{image.default}}
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}. 


}

\examples{
data("dna.obj")
x<-dna.obj

## Display only polymorphic sites without gaps
image(x,all=FALSE,fifth=FALSE,show.names=TRUE,cex=0.6)

## Display only polymorphic sites with gaps
image(x,all=FALSE,fifth=TRUE,show.names=TRUE,cex=0.6)

\dontrun{
## Display entire sequences
image(x,all=FALSE,show.names=TRUE,cex=0.6)
}
}


\keyword{DNA}

