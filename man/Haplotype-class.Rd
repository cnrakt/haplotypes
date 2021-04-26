\name{Haplotype-class}
\Rdversion{1.1}
\docType{class}
\alias{Haplotype}
\alias{Haplotype-class}
\title{Class \code{"Haplotype"} in the Package \pkg{haplotypes}}
\description{
S4 class to store haplotype information.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Haplotype", haplist, hapind, uniquehapind, sequence, d, freq, nhap)},  however use function \code{\link{haplotype}}  instead.
}
\section{Slots}{
  \describe{
    \item{\code{haplist}:}{Object of class \code{"list"}, containing the names of individuals that share the same haplotype. }
    \item{\code{hapind}:}{Object of class \code{"list"}, containing the index of individuals that share the same haplotype. }
    \item{\code{uniquehapind}:}{Object of class \code{"numeric"}, containing the index of the first occurrence of unique haplotypes. }
    \item{\code{sequence}:}{Object of class \code{"matrix"} if present, giving the DNA sequence matrix of unique haplotypes. }
    \item{\code{d}:}{Object of class \code{"matrix"}, giving the absolute pairwise character difference matrix of unique haplotypes. }
    \item{\code{freq}:}{Object of class \code{"numeric"}, giving the haplotype frequencies. }
    \item{\code{nhap}:}{Object of class \code{"numeric"}, giving the total number of haplotypes.}
  }
}
\section{Methods}{
  \describe{
    \item{as.dna}{\code{signature(x = "Haplotype")}: if Haplotype object contains dna sequences, coerces an object of class Haplotype to an object of class Dna, else returns an error message.}
    \item{as.list}{\code{signature(x = "Haplotype")}: assigns slots of an object Haplotype to list elements.}
    \item{grouping}{\code{signature(x = "Haplotype")}: creates a matrix with haplotypes as rows, grouping factor (populations, species, etc.) as columns and abundance as entries. }
    \item{hapreord}{\code{signature(x = "Haplotype")}: reorders haplotypes according to the ordering factor.}
    \item{length}{\code{signature(x = "Haplotype")}: returns the number of haplotypes. }
     \item{pieplot}{\code{signature(x = "Parsimnet", y = "Haplotype")}: plot pie charts on statistical parsimony network.}
     \item{pielegend}{\code{signature(x = "Parsimnet", y = "Haplotype")}: add legends to pie charts produced using \code{\link{pieplot}}.}
    \item{show}{\code{signature(object = "Haplotype")}: displays the object briefly. }
	 }
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}


\keyword{HAPLOTYPE ANALYSIS}
\keyword{CLASSES}
