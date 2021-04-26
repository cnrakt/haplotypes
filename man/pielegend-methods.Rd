\name{pielegend-methods}
\docType{methods}
\alias{pielegend}
\alias{pielegend-methods}
\alias{pielegend,Parsimnet,Haplotype-method}
\title{Add Legends to Plots}

\description{
	This function can be used to add legends to pie charts produced using \code{\link{pieplot}}.
}

\usage{
\S4method{pielegend}{Parsimnet,Haplotype}(p,h,net=1,factors,...)

}

\arguments{
  \item{p}{ a \code{\link{Parsimnet}} object. }
  \item{h}{ a \code{\link{Haplotype}} object.}
  \item{net}{a numeric vector of length one indicating which network to plot.}
  \item{factors}{a vector or factor giving the grouping variable (populations, species, etc.), with one element per individual. }
  \item{...}{ arguments to be passed to \code{\link{legend}} and others. See \sQuote{Details} }

}


\details{
This method calls \code{\link{legend}} \code{{graphics}}, some default parameters changed:

  \describe{
  	
  	\item{col}{ an integer or character vector for the edge colors. By default, it is \code{rainbow}. }
  	 \item{fill}{ an integer or character vector for the filling colors. By default, it is \code{rainbow} }
	\item{legend}{ a numeric or character vector to appear in the legend. By default, it is the levels of the grouping factor for haplotypes. }
  	\item{x}{ position of the legend. Default is set to \code{"topright"}. }
  	 

  }
}

 \value{
 See \sQuote{ \code{\link{legend}} \code{{graphics}}}
}

\section{Methods}{
\describe{

\item{\code{signature(p = "Parsimnet", h = "Haplotype")}}{
}
}}

\seealso{
 \code{\link{plot,Parsimnet-method}},\code{\link{floating.pie}}, \code{\link{plot.default}}, \code{\link{plot.network.default}} and \code{\link{legend}}.
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}. 


}


\examples{
data("dna.obj")
x<-dna.obj
h<-haplotypes::haplotype(x)

### Statistical parsimony with 95% connection limit 
p<-parsimnet(x) 

#randomly generated populations
pop<-c("pop1","pop2","pop3","pop4","pop5","pop6","pop7","pop8")
set.seed(5)
pops<-sample(pop,nrow(x),replace=TRUE)

## Plotting with default parameters.
pieplot(p,h,1, pops)

## Add legend with default parameters.
pielegend(p,h,1,pops)

## Change colors for the populations.
#8 colors for 8 populations 
cols<-colors()[c(30,369,552,558,538,642,142,91)]
pieplot(p,h,1, pops,col=cols)
pielegend(p,h,1,pops,col= cols)

}

\keyword{STATISTICAL PARSIMONY}