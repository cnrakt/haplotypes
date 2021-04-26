\name{Parsimnet-class}
\Rdversion{1.1}
\docType{class}
\alias{Parsimnet}
\alias{Parsimnet-class}

\title{Class \code{"Parsimnet"} in the Package \pkg{haplotypes} }
\description{
S4 class to store statistical parsimony networks and additional information.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Parsimnet", d, tempProbs, conlimit, prob, nhap, rowindex), however use function \link{parsimnet} instead}.

}
\section{Slots}{
  \describe{
    \item{\code{d}:}{Object of class \code{"list"} containing the geodesic distance matrix of haplotypes and intermediates for each network. }
    \item{\code{tempProbs}:}{Object of class \code{"numeric"} giving the probabilities of parsimony for mutational steps beyond the connection limit.  }
    \item{\code{conlimit}:}{Object of class \code{"numeric"} giving the number of maximum connection steps at connection limit.}
    \item{\code{prob}:}{Object of class \code{"numeric"} giving the user defined connection limit.  }
    \item{\code{nhap}:}{Object of class \code{"numeric"} giving the number of haplotypes in each network.}
    \item{\code{rowindex}:}{Object of class \code{"list"} containing vectors giving the index of haplotypes in each network. }
  }
}

\section{Methods}{
  \describe{
    \item{as.list}{\code{signature(x = "Parsimnet")}: assigns slots of an object Parsimnet to list elements.}
    \item{as.network}{\code{signature(x = "Parsimnet")}: coerces \code{Parsimnet} object to \code{\link{network}} \code{{network}} object }
    \item{as.networx}{\code{signature(x = "Parsimnet")}: coerces \code{Parsimnet} object to \code{\link{networx}} \code{{phangorn}} object  }
    \item{length}{\code{signature(x = "Parsimnet")}: returns the length of network(s). }
    \item{names}{\code{signature(x = "Parsimnet")}: gets names of networks in \code{Parsimnet} object }
    \item{names<-}{\code{signature(x = "Parsimnet")}: sets names of networks in \code{Parsimnet} object }
    \item{plot}{\code{signature(x = "Parsimnet")}: plots statistical parsimony networks.}
    \item{pieplot}{\code{signature(x = "Parsimnet", y = "Haplotype")}: plots pie charts on statistical parsimony networks }
    \item{pielegend}{\code{signature(x = "Parsimnet", y = "Haplotype")}: add legends to pie charts produced using \code{\link{pieplot}}.}
    \item{rownames}{\code{signature(x = "Parsimnet")}: gets names of vertices in networks. }
    \item{rownames<-}{\code{signature(x = "Parsimnet")}: gets names of vertices in networks }
    \item{show}{\code{signature(object = "Parsimnet")}: displays the object briefly. }
	 }
}




\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\keyword{STATISTICAL PARSIMONY}
\keyword{CLASSES}

