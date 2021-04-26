\name{haplotypes-package}
\alias{haplotypes-package}
\alias{haplotypes}
\docType{package}
\title{
Manipulating DNA Sequences and Estimating Unambiguous Haplotype Network with Statistical Parsimony
}
\description{
This package provides S4 classes and methods for reading and manipulating aligned DNA sequences, supporting an indel coding methods (only simple indel coding method is available in the current version), showing base substitutions and indels, calculating absolute pairwise distances between DNA sequences, and collapses identical DNA sequences into haplotypes or infering haplotypes using user provided absolute pairwise character difference matrix.  This package also includes S4 classes and methods for estimating genealogical relationships among haplotypes using statistical parsimony and plotting parsimony networks. 
}

\author{
Caner Aktas, \email{caktas.aca@gmail.com}
}

\examples{

## Read example FASTA file.	
f<-system.file("example.fas",package="haplotypes") 
# invalid character 'N' was replaced with '?' with a warning message
x<-read.fas(file=f)
# an object of class 'Dna'
x 

## or load DNA Sequence data set.
data("dna.obj") 
x<-dna.obj
\dontrun{
x
}

## Compute an absolute pairwise character difference matrix from DNA sequences.
# coding gaps using simple indel coding method
d<- distance(x,indels="sic") 
\dontrun{
d
}

## Infer haplotypes using the 'Dna' object.
# coding gaps using simple indel coding method
h<-haplotype(x,indels="s") 
\dontrun{
h
}

## Conduct statistical parsimony analysis with 95% connection limit.
#algortihmic method
\dontrun{
p<-parsimnet(x,prob=.95) 
p
# plot network
plot(p) 
}

## Plotting pie charts on the statistical parsimony network
\dontrun{
data("dna.obj")
x<-dna.obj
h<-haplotypes::haplotype(x)

## Statistical parsimony with 95% connection limit 
p<-parsimnet(x) 

#randomly generated populations
pop<-c("pop1","pop2","pop3","pop4","pop5","pop6","pop7","pop8")
set.seed(5)
pops<-sample(pop,nrow(x),replace=TRUE)

# Plotting with default parameters.
pieplot(p,h,1, pops)

}
}

\keyword{PACKAGE}
