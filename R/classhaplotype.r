###Class Haplotype, functions and methods
##march 2015, cnrakt


#Class Haplotype
setClass(Class="Haplotype", representation=representation(haplist ="list",hapind ="list", uniquehapind="numeric",sequence="matrix",d="matrix",freq="numeric",nhap="numeric"))


#Set initialize Class Haplotype
setMethod("initialize", "Haplotype", function(.Object,haplist= list(),hapind =list(),uniquehapind=numeric(),sequence=matrix(,0,0),d=matrix(,0,0),freq=sapply(.Object@haplist,length),nhap=length(.Object@haplist)) {
	
	.Object@haplist<-haplist
	.Object@hapind<-hapind
	.Object@uniquehapind<-uniquehapind
	.Object@sequence<-sequence
	.Object@d<-d
	.Object@freq<-freq
	.Object@nhap<-nhap
	.Object
})


#Show method for Haplotype objects

setMethod("show","Haplotype", function(object)
{
	
	cat("*** S4 Object of Class Haplotype ***\n\n")
	
	show(object@haplist)
	cat("\nNumber of haplotypes : ", object@nhap)
	cat("\n\nslots of an object Haplotype:\n")
	cat(slotNames(object),"\n")
	
})


#internal function: inferring haplotypes from distance matrix or dist object

haplo<-function(x)
{
	
	x <- as.matrix(x)
	diag(x)<-0
	nseq<-nrow(x)
	whap<-x[1,]==0
	haploind<-list(which(whap))
	haplovec<-which(whap)
	
	if(length(x)>1)
	{ 
		for(i in 2:nseq)
		{ 
			whap<-x[i,]==0
			whap[haplovec]<-FALSE
			haploind<-c(haploind,list(which(whap)))
			haplovec<-unique(c(haplovec,which(whap)))
		}	
	}
	empthaplo<-sapply(haploind,length)
	haploind<-unique(haploind[empthaplo>0])
	haplolist<-lapply(haploind,names)
	
	uniqhapindex<-sapply(haploind,"[",1)
	
	hapdistmat<-as.matrix(x[uniqhapindex,uniqhapindex])
	
	
	freq<-sapply(haplolist,length)
	hapnum<-length(freq)
	
	names(haploind)<-paste("haplotype",1:hapnum,sep="")
	names(haplolist)<-names(haploind)
	
	
	hapobj<-new("Haplotype",haplist=haplolist,hapind=haploind,uniquehapind=uniqhapindex,d=hapdistmat,freq=freq,nhap=hapnum)
	
	hapobj
	
}



#Generic haplotype
setGeneric (
name= "haplotype",
def=function(x,...)standardGeneric("haplotype")
)


#Dna object to haplotype object

setMethod("haplotype","Dna", function(x,indels="sic")
{
	
	d<-distance(x,indels=indels)
	if(length(d)==0) stop("at least two DNA sequences are required")
	h<-haplo(d)
	h@sequence<-x@sequence[h@uniquehapind,,drop=FALSE]
	h
})


#matrix object to haplotype object

setMethod("haplotype","matrix", function(x)
{
	if(nrow(x)==1) stop("dimension of the distance matrix (x) must be greater than one")
	haplo(x)	
})


#dist object to haplotype object

setMethod("haplotype","dist", function(x)
{
	if(length(x)==0) stop("length of the distance object (x) must be greater than zero")
	haplo(x)	
})



#Coerce Haplotype object to a list

setMethod("as.list","Haplotype", function(x) 
{
	l<-list(x@haplist,x@hapind, x@uniquehapind,x@sequence,x@d,x@freq,x@nhap)
	names(l)<-c("haplist","hapind","uniquehapind","sequence","d","freq","nhap")
	l	
})



#length method for Dna objects

setMethod("length","Haplotype", function(x) 
{
	x@nhap		
})


#Generic hapreord

setGeneric (
name= "hapreord",
def=function(x,...)standardGeneric("hapreord")
)


#hapreord method for Haplotype objects

setMethod(f="hapreord", signature= "Haplotype", definition=function(x,order=c(1:x@nhap))
{
	if(length(order)!=length(x))  stop(paste("'order' must be vector of length",length(x)))
	if(any(order<1)|any(order>length(x))) stop(paste("elements of 'order' must be integers in the range [1,",length(x),"]",sep=""))
	haplolist<-x@haplist[order]
	names(haplolist)<-paste("haplotype",1:x@nhap,sep="")
	haploind<-x@hapind[order]
	names(haploind)<-names(haplolist)
	
	if(nrow(x@sequence)>0) 
	{
		hapobj<-new("Haplotype",haplist=haplolist,hapind=haploind,uniquehapind=x@uniquehapind[order],sequence=x@sequence[order,,drop=FALSE],d=x@d[order,order],freq=x@freq[order])
		}	else {
			hapobj<-new("Haplotype",haplist=haplolist,hapind=haploind,uniquehapind=x@uniquehapind[order],d=x@d[order,order],freq=x@freq[order])
			
		} 
	
	return(hapobj)
})



#Argument 'factor' is renamed, Dec 2017
#Generic grouping

setGeneric (
name= "grouping",
def=function(x,...)standardGeneric("grouping")
)


#grouping method for Haplotype objects

setMethod(f="grouping", signature= "Haplotype", definition=function(x,factors)
{
	l<-x@nhap
	flevels<-levels(factor(factors))
	hapmat<-matrix(0, l, length(flevels))
	hapvec<-vector("numeric",length(factors))
	rownames(hapmat)<-1:l
	colnames(hapmat)<-flevels
	
	for(i in 1:l)
	{
		hind<- x@hapind[[i]]
		fstab<-table(factor(factors[hind]))
		hapmat[i,names(fstab)]<-fstab
		hapvec[hind]<-i
	}
	return(list(hapmat=hapmat,hapvec=hapvec))
	
})





#fixing a bug in seqlengths, May 2016
#Coerce Haplotype objects to Dna object

setMethod(f="as.dna", signature= "Haplotype", definition=function(x)
{
    
    seq<-x@sequence
    if(!nrow(seq)) stop("Haplotype object does not contain DNA sequences")
    dnaobj<-new("Dna",sequence=seq,seqlengths=rep(ncol(seq),nrow(seq)),seqnames=rownames(seq))
    dnaobj
}
)



















