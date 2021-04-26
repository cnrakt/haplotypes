\name{read.fas}
\alias{read.fas}
\title{
Read sequences from a file in FASTA format
}
\description{
Read DNA sequences from a file in FASTA Format.
}
\usage{
read.fas(file)
}
\arguments{
  \item{file}{
  	the name of the file, which the sequence in the FASTA format is to be read from. If it does not contain an \emph{absolute} path, the file name is \emph{relative} to the current working directory,
    \code{\link{getwd}()}.
}
}

\value{
\code{read.fas} returns an object of class  \code{\link{Dna}}.

}


\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}
\note{
By default, valid characters are "A","C","G","T","a","c","g","t","-","?" for the class  \code{\link{Dna}}. Numeric entries (integers) between 0-5 will be converted to "?","A","C","G","T","-", respectively. Invalid characters will be replaced with "?" with a warning message. 
}

\seealso{
\code{\link{Dna}}
}
\examples{

##Reading example file.	
f<-system.file("example.fas",package="haplotypes")
 
# invalid character 'N' was replaced with '?' with a warning message
x<-read.fas(file=f) 

# an object of class 'Dna'
x

}


\keyword{DNA ANALYSIS}
