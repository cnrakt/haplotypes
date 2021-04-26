\name{as.dna-methods}
\docType{methods}
\alias{as.dna}
\alias{as.dna-methods}
\alias{as.dna,character-method}
\alias{as.dna,data.frame-method}
\alias{as.dna,DNAbin-method}
\alias{as.dna,Haplotype-method}
\alias{as.dna,list-method}
\alias{as.dna,matrix-method}
\alias{as.dna,phyDat-method}



\title{Coerces an object to a \code{Dna} object}
\description{
 Coerces an object that contains DNA sequences to an object of Class \code{\link{Dna}}.
 
}
\usage{
\S4method{as.dna}{matrix}(x)
\S4method{as.dna}{data.frame}(x)
\S4method{as.dna}{list}(x)
\S4method{as.dna}{character}(x)
\S4method{as.dna}{Haplotype}(x)
\S4method{as.dna}{DNAbin}(x)
\S4method{as.dna}{phyDat}(x)

}

\arguments{
  \item{x}{ a matrix, a data.frame, a list, a character, an object of class \code{\link{Haplotype}}, \code{\link{DNAbin}} \code{{ape}} or \code{\link{phyDat}} \code{{phangorn}} object containing the DNA sequences.}
}

\details{
	
Elements of the list must be vectors. Each element of the list contains a single DNA sequence. If the sequence lengths differ, the longest sequence is taken into account and gaps are introduced to the shorter sequences at the end of the matrix in the slot \code{sequence}. Sequence length information is stored in the slot \code{seqlengths}.

Valid characters for the slot \code{sequence} are "A","C","G","T","a","c","g","t","-","?". During the conversion of the object to the class \code{\link{Dna}}, integers 0,1,2,3,4,5 or characters "0","1","2","3","4","5" are converted to "?","A","C","G","T","-", respectively. Invalid characters are replaced with "?" with a warning message.  
  
}


\value{
  an object of class \code{\link{Dna}.}
  }


\section{Methods}{
\describe{

\item{\code{signature(x = "matrix")}}{
coerces matrix to a \code{Dna} object.
}

\item{\code{signature(x = "data.frame")}}{
coerces data.frame to a \code{Dna} object.
}

\item{\code{signature(x = "list")}}{
coerces list to a \code{Dna} object.
}

\item{\code{signature(x = "character")}}{
coerces characters to a \code{Dna} object.
}

\item{\code{signature(x = "Haplotype")}}{
coerces a \code{\link{Haplotype}} object to a \code{Dna} object.
}

\item{\code{signature(x = "DNAbin")}}{
coerces a \code{\link{DNAbin}} object to a \code{Dna} object.
}

\item{\code{signature(x = "phyDat")}}{
coerces a \code{\link{phyDat}} object to a \code{Dna} object.
}

}}


\examples{

## Coercing a matrix to a 'Dna' object.
# all valid characters
x<-matrix(c("?","A","C","g","t","-","0","1","2","3","4","5"),4,6)
rownames(x)<-c("seq1","seq2","seq3","seq4") 
dna.obj<-as.dna(x) 
dna.obj
# the sequence matrix
dna.obj@sequence

\dontrun{
# includes invalid characters
x<-matrix(c("X","y","*","?","t","-","0","1","2","3","4","5"),4,6) 
rownames(x)<-c("seq1","seq2","seq3","seq4") 
dna.obj<-as.dna(x) 
dna.obj
dna.obj@sequence

# all valid integers
x<-matrix(c(0,1,2,3,4,5,0,1,2,3,4,5),4,6) 
rownames(x)<-c("seq1","seq2","seq3","seq4") 
dna.obj<-as.dna(x) 
dna.obj
dna.obj@sequence

## Coercing a data.frame to a 'Dna' object.
x<-data.frame(matrix(c("?","A","C","g","t","-","0","1","2","3","4","5"),4,6)) 
rownames(x)<-c("seq1","seq2","seq3","seq4") 
dna.obj<-as.dna(x) 
dna.obj
dna.obj@sequence

## Coercing a list to a 'Dna' object.
seq1<-c("?","A","C","g","t","-","0","1")
seq2<-c("?","A","C","g","t","-","0","1","2")
seq3<-c("?","A","C","g","t","-","0","1","2","3")
x<-list(seq1=seq1,seq2=seq2,seq3=seq3) 
dna.obj<-as.dna(x)

# sequence lengths differ 
dna.obj@seqlengths
dna.obj@sequence

## Coercing a character vector to a Dna object.
seq<-c("?","A","C","g","t","-","0","1")
x<-as.dna(seq)
x


## Coercing a Haplotype object to a Dna object.
data("dna.obj")
x<-dna.obj
h<-haplotype(x)

# DNA Sequences of unique haplotypes  
dna.obj<-as.dna(h)
dna.obj

d<-distance(x)

# if 'Haplotype' object does not contain 'DNA' Sequences
h<-haplotype(d)

# returns an error
as.dna(h)

## Coercing a DNAbin object to a Dna object.
require(ape)
data(woodmouse)
x<-as.dna(woodmouse)
x

}}

\keyword{DNA ANALYSIS}
