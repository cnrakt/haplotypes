\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,Dna-method}
\alias{show,Haplotype-method}
\alias{show,Parsimnet-method}
\title{Methods for function \code{show} in the package \pkg{haplotypes}}
\description{
Show objects of classes \code{Dna}, \code{Haplotype}, and  \code{Parsimnet} 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "Dna")}}{
displays \code{Dna} object briefly: The total number of DNA sequences, names of the first six sequences (if nrow(x)>=6), length of the shortest and longest sequences and the names of the slots.
}

\item{\code{signature(object = "Haplotype")}}{
displays \code{Haplotype} object briefly: The list of individuals that share the same haplotypes, the total number of haplotypes and the names of the slots.
}

\item{\code{signature(object = "Parsimnet")}}{
displays \code{Parsimnet} object briefly: The total number of networks, the maximum connection steps at chosen probability, the total number of haplotypes in each network,the total number of intermediates in each network, total network lengths (scores) of each network and the names of the slots.}

}}

\keyword{DNA ANALYSIS}
\keyword{HAPLOTYPE ANALYSIS}
\keyword{STATISTICAL PARSIMONY}
