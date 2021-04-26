\name{pieplot-methods}
\docType{methods}
\alias{pieplot}
\alias{pieplot-methods}
\alias{pieplot,Parsimnet,Haplotype-method}
\title{Plots pie charts on statistical parsimony network}
\description{
	Plotting pie charts on the statistical parsimony network.
}

\usage{
\S4method{pieplot}{Parsimnet,Haplotype}(x,y,net=1,factors, coord = NULL,inter.labels=FALSE,interactive=FALSE,rex=1,...)

}

\arguments{
  \item{x}{ a \code{\link{Parsimnet}} object. }
  \item{y}{ a \code{\link{Haplotype}} object.}
  \item{net}{a numeric vector of length one indicating which network to plot.}
  \item{factors}{a vector or factor giving the grouping variable (populations, species, etc.), with one element per individual. }
  \item{coord}{a matrix that contains user specified coordinates of the vertices, or NULL to  generate vertex layouts with \code{\link{network.layout}} function from package \pkg{network}. }
  \item{inter.labels}{boolean; should vertex labels of intermediate haplotypes be displayed? }
  \item{interactive}{boolean; should vertices be interactively adjusted? }
  \item{rex}{ expansion factor for the pie radius.}
  \item{...}{ arguments to be passed to \code{\link{floating.pie}} and others. See \sQuote{Details} }
}

\details{
This method calls \code{\link{floating.pie}} \code{{plotrix}}, \code{\link{network.vertex}} \code{{network}}, and \code{\link{plot.default}}, \code{\link{lines}}, and \code{\link{text}}  \code{{graphics}}. This method also uses some internal structures of \code{\link{plot.network.default}} from package \pkg{network}. The following additional arguments can be passed to these functions: 

  \describe{
  	 \item{mode}{ the vertex placement algorithm. Default is set to \code{"fruchtermanreingold"}. }
     \item{pad}{ amount to pad the plotting range; useful if labels are being clipped. Default is set to 1.}
     \item{displaylabels}{boolean; should vertex labels be displayed? }
   \item{label}{ a vector of vertex labels. By default, the rownames of the distance matrix (\code{rownames(p@d[[net]])}) are used. If \code{inter.labels==FALSE} only haplotype labels are displayed. }
       \item{label.cex}{ character expansion factor for labels. Default is set to 0.75.}
       \item{label.col}{ an integer or character vector for the label colors. By default, it is 1 (black). }
       \item{label.pos}{ position at which labels should be placed relative to vertices. 0 and 6 results in labels which are placed away from the center of the plotting region; 1, 2, 3, and 4 result in labels being placed below, to the left of, above, and to the right of vertices, respectively; and label.pos 5 or greater than 6 results in labels which are plotted with no offset (i.e., at the vertex positions). Default is set to 0.}
       \item{label.pad}{ amount to pad the labels. This setting is available only if the labels are plotted with offset relative to vertex positions. Default is set to 1.}
       \item{vertex.cex}{ a numeric vector of expansion factor for intermediate vertices (only). By default it is \code{(0.5)*min(radius)}. Use 'radius' to specify size of pie charts. }
  	\item{col}{ the colors of the pie sectors (i.e., colors for populations), by default \code{rainbow}. }
  	\item{vertex.col}{ an integer or character vector for the intermediate vertex colors. By default, it is 1 (black). }
	\item{edge.col}{ an integer or character vector for the edge colors. By default, it is 1 (black). }
	 	\item{edge.lwd}{ a numeric vector, edges line width. By default, it is 1. }
	 	\item{edge.lty}{ a numeric vector of length one, specifies the line type for the edges. By default it is 1. }
		\item{edges}{ the number of lines forming a pie circle, By default, it is 200. }
		\item{radius}{ a numeric vector of length \code{p@nhap[net]} for the radius of drawn pie charts. Useful for specifying the radius independent of the haplotype frequencies. Default is \code{(0.8*(haplotype frequencies)*rex)/max(haplotype frequencies)}. }
		\item{vertex.sides}{ number of polygon sides for vertices. Default is set to 50. }
	\item{xlab}{x axis label. }
	\item{xlab}{y axis label. }

  }
}

 \value{
 A two-column matrix containing the vertex positions as x,y coordinates.
}


\section{Methods}{
\describe{

\item{\code{signature(x = "Parsimnet", y = "Haplotype")}}{
}
}}


  \note{ 
  Some internal structures of \code{\link{plot.network.default}} is taken from package \pkg{network} with modifications, author Carter T. Butts.
}

\seealso{
 \code{\link{plot,Parsimnet-method}},\code{\link{floating.pie}}, \code{\link{plot.default}} and \code{\link{plot.network.default}}
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

## Change colors for the populations.
#8 colors for 8 populations 
cols<-colors()[c(30,369,552,558,538,642,142,91)]
pieplot(p,h,1, pops,col=cols)

## Expanding pie charts and intermediate vertices.
pieplot(p,h,1, pops,rex=2)

## Adjusting intermediate vertex sizes.
pieplot(p,h,1, pops, vertex.cex=rep(0.2, nrow(p@d[[1]])-p@nhap))

## Expanding pie charts and intermediate vertices, adjusting intermediate vertex sizes.
pieplot(p,h,1, pops,rex=2, vertex.cex=rep(0.1, nrow(p@d[[1]])-p@nhap))

## Adjusting radius of pie charts.
pieplot(p,h,1, pops,radius=rep(1, p@nhap))

\dontrun{ 
## Interactively adjusting vertex positions.
pieplot(p,h,1, pops, interactive=TRUE)
}

### Multiple networks with 99% connection limit.
p<-parsimnet(x,prob=.99) 

## Plotting first network with default parameters.
pieplot(p,h,1, pops)

## Change colors for the populations.
#8 colors for 8 populations
cols<-colors()[c(30,369,552,558,538,642,142,91)]
pieplot(p,h,1, pops,col=cols)







}

\keyword{STATISTICAL PARSIMONY}

