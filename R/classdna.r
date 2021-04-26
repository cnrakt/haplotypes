### class Dna, functions and methods
##march 2015, cnrakt

##april 2016, cnrakt
##An internal code Statphi is taken from package ade4 version 1.7-8 without any modification, author Sandrine Pavoine


#class Dna
setClass(Class="Dna", representation=representation(sequence="matrix",seqlengths="numeric",seqnames="character"))


#validity: testing Dna object

validity <- function(x)
{
	seq<-x@sequence
	

	seq[seq==0]<-"?"
	seq[seq==1]<-"A"
	seq[seq==2]<-"C"
	seq[seq==3]<-"G"
	seq[seq==4]<-"T"
	seq[seq==5]<-"-"
	
	unchar<- setdiff(c(seq),c("A","C","G","T","a","c","g","t","-","?"))
	
	if(length(unchar)>0)
	{
		
		seq[!is.na(match(c(seq),unchar))]<-"?"
		if(length(unchar)==1) warning(paste("invalid character ","'", paste(unchar,collapse=","),"'", " was replaced with '?'",sep=""),"\n(valid characters: ", paste(c("A","C","G","T","a","c","g","t","-","?","0","1","2","3","4","5"),collapse=","),")")
		if(length(unchar)>1) warning(paste("invalid characters ","'", paste(unchar,collapse=","),"'", " were replaced with '?'",sep=""),"\n(valid characters: ", paste(c("A","C","G","T","a","c","g","t","-","?","0","1","2","3","4","5"),collapse=","),")")
	}

	x@sequence<-seq
	x
	
}


#Set initialize Class Dna
setMethod("initialize", "Dna", function(.Object,sequence=matrix(,0,0),seqlengths=integer(0),seqnames=character(0)) {
	
	.Object@sequence<-sequence
	.Object@seqlengths<-seqlengths
	if(!is.character(seqnames)) seqnames<-as.character(seqnames)
	if(nrow(sequence)>0&length(seqnames)==0) seqnames<-as.character(1:nrow(sequence))
	if(nrow(sequence)>0&length(seqlengths)==0) .Object@seqlengths<-rep(ncol(sequence),nrow(sequence))
	.Object@seqnames<-seqnames
	validity(.Object)
})




#Show method for Dna objects
setMethod("show","Dna", function(object)
{
	
	cat("*** S4 Object of Class Dna ***\n\n")
	cat("\nNumber of DNA sequences: \n")
	cat(nrow(object@sequence),"\n")
	cat("\nNames: \n")
	if(nrow(object@sequence)>6) cat(c(head(object@seqnames),"..."),"\n") else cat(head(object@seqnames),"\n")
	cat("\nLength of Shortest DNA Sequence: \n")
	if(length(object@seqlengths)==0) cat(integer(0)) else cat(min(object@seqlengths),"\n")
	cat("\n\nLength of Longest DNA Sequence: \n")
	if(length(object@seqlengths)==0) cat(integer(0)) else cat(max(object@seqlengths),"\n")
	cat("\n\nslots of an object Dna:\n")
	cat(slotNames(object),"\n")
	
})


#Coerce Dna objects to matrix
setMethod("as.matrix","Dna", function(x)
{
	x@sequence	
})


#Coerce Dna objects to data.frame
setMethod("as.data.frame","Dna", function(x)
{
	df<-as.data.frame(x@sequence)
	colnames(df)<-1:ncol(df)
	df
	
})


#Coerce Dna objects to list
setMethod("as.list","Dna", function(x) 
{
	n<-nrow(x@sequence)
	l<-vector("list",n)
	lseq<-x@seqlengths
	for(i in 1:n) l[[i]]<-x@sequence[i,1:lseq[i]]
	names(l)<-x@seqnames
	l	
})

#Coerce Dna objects to numeric matrix

setMethod("as.numeric","Dna", function(x)
{
	seq<-toupper(x@sequence)
	seq[seq=="?"]<-0
	seq[seq=="A"]<-1
	seq[seq=="C"]<-2
	seq[seq=="G"]<-3
	seq[seq=="T"]<-4
	seq[seq=="-"]<-5
	x<-matrix(as.numeric(unlist(seq)),nrow=nrow(seq),ncol=ncol(seq),dimnames=list(x@seqnames,1:ncol(seq)))
	x
})	


#names method for Dna objects

setMethod("names","Dna", function(x) 
{
	x@seqnames		
})


#names replace method for Dna objects

setReplaceMethod("names","Dna", function(x,value) 
{
	if(is.numeric(value)) value<-as.character(value)
	rownames(x@sequence)<-value
	x@seqnames<-value
	x
})


#length method for Dna objects

setMethod("length","Dna", function(x) 
{
	ncol(x@sequence)		
})


#ncol method for Dna objects

setMethod("ncol","Dna", function(x) 
{
	ncol(x@sequence)		
})



#nrow method for Dna objects

setMethod("nrow","Dna", function(x) 
{
	nrow(x@sequence)		
})



#fixing a bug in Extract, July 2016
#Extract method for Dna objects

setMethod("[","Dna", function(x,i=1:nrow(x),j=1:ncol(x),as.matrix=TRUE)
{
    seq<-x@sequence[i,j,drop=FALSE]
    lseq<-ncol(seq)
    seqnames<-x@seqnames[i]
    if(as.matrix) return(seq) else new("Dna",sequence=seq,seqlengths=rep(lseq,nrow(seq)),seqnames=seqnames)
    
})


#Extract replace method for Dna objects

setReplaceMethod("[","Dna", function(x,i,j,value) 
{
	x@sequence[i,j]<-value
	x<-validity(x)
	x
})


#range method for Dna objects
setMethod("range","Dna", function(x)
{
    range(x@seqlengths)
})


#tolower method for Dna objects
setMethod("tolower","Dna", function(x)
{
   
    x@sequence<-tolower(x@sequence)
 x
})
#toupper method for Dna objects
setMethod("toupper","Dna", function(x)
{
    
    x@sequence<-toupper(x@sequence)
	x
})

#unique method for Dna objects
setMethod("unique","Dna", function(x,gaps=FALSE)
{
    x<-toupper(x)
    xl<-as.list(x)
    
    if(!gaps)
    {
        for(i in 1:nrow(x))
        {
            r<-xl[[i]]
            xl[[i]]<-r[r!="-"]
        }
    }
    unique(xl)
    
})


#rownames method for Dna objects

setMethod("rownames","Dna", function(x)
{
    rownames(x@sequence)
})


#rownames replace method for Dna objects

setReplaceMethod("rownames","Dna", function(x,value)
{
    if(is.numeric(value)) value<-as.character(value)
    rownames(x@sequence)<-value
    x@seqnames<-value
    x
})





#Generic as.dna

setGeneric (
name= "as.dna",
def=function(x,...)standardGeneric("as.dna")
)


#Coerce matrix to Dna object

setMethod(f="as.dna", signature= "matrix", definition=function(x)
{
	if(is.numeric(x))
	{
		x[x==0]<-"?"
		x[x==1]<-"A"
		x[x==2]<-"C"
		x[x==3]<-"G"
		x[x==4]<-"T"
		x[x==5]<-"-"
		x[x==NA]<-NA
	}
	
	dnaobj<-new("Dna",sequence=x,seqlengths=rep(ncol(x),nrow(x)),seqnames=rownames(x))
	dnaobj
}
)



#Coerce data.frame to Dna object

setMethod(f="as.dna", signature= "data.frame", definition=function(x)
{
	if(is.numeric(x))
	{
		x[x==0]<-"?"
		x[x==1]<-"A"
		x[x==2]<-"C"
		x[x==3]<-"G"
		x[x==4]<-"T"
		x[x==5]<-"-"
		x[x==NA]<-NA
	}
	
	dnaobj<-new("Dna",sequence=as.matrix(x),seqlengths=rep(ncol(x),nrow(x)),seqnames=rownames(x))
	dnaobj
}
)



#Coerce list to Dna object

setMethod(f="as.dna", signature= "list", definition=function(x)
{

	lseq<-sapply(x,length)
	l<-length(x)
	seqnames<-names(x)
	maxseql<-max(lseq)
	
	seq.mat<-matrix("-",l,maxseql,dimnames=list(seqnames,1:maxseql))
	for(i in 1:(l)) seq.mat[i,1:lseq[i]]<-x[[i]]
	
	dnaobj<-new("Dna",sequence=seq.mat,seqlengths=lseq,seqnames=seqnames)
	dnaobj
}
)


#Coerce character vector to Dna object
setMethod(f="as.dna", signature= "character", definition=function(x)
{
    lseq<-length(x)
    seqnames<-names(x)
    
    seq.mat<-matrix("-", 1, lseq,dimnames=list(seqnames,1: lseq))
    seq.mat[1,1:lseq]<-x
    
    dnaobj<-new("Dna",sequence=seq.mat,seqlengths=lseq,seqnames=seqnames)
    dnaobj
}
)
#setOldClass("DNAbin")
#Coerce DNAbin objects to Dna object
setMethod(f="as.dna", signature= "DNAbin", definition=function(x)
{
    
    as.dna(as.character(x))
}
)


#Coerce phyDat objects to Dna object
setMethod(f="as.dna", signature= "phyDat", definition=function(x)
{
    
    as.dna(as.character(x))
}
)




#Bugs corrected, gives less errors, Dec 2017 
#internal function: read fasta, the first line of the files must start with a ">" (greater-than) symbol.

read.fas<-function(file)

{
	fas<-scan(file, what="char",quote="", sep="\n", strip.white=TRUE, quiet= TRUE)
	
	l<-length(fas)
	if(!l) stop(paste(file,"file is empty."))
	heads<-grep(">",fas)
	if(!length(heads)) stop("invalid fasta format, description lines of each sequence must start with a '>' (greater-than) symbol")
	tes<- diff(heads)
	test1<-any(tes==1)
	
	if(test1) stop("invalid fasta format, description lines of each sequence must start with a '>' (greater-than) symbol and file must not contain empty sequences")

	seqnames<-gsub(">","",fas[heads])
	len<-length(seqnames)
	seqlist<-vector("list", len)
	names(seqlist)<-seqnames
	for(i in 1:(len-1))
	{
		seqlist[[i]]<-strsplit(paste(fas[(heads[i]+1):(heads[i+1]-1)],collapse=""),split="")[[1]]
		
	}
	seqlist[[i+1]]<-strsplit(paste(fas[(heads[i+1]+1):(l)],collapse=""),split="")[[1]]
	
	as.dna(seqlist)
		

}




#Generic subs

setGeneric (
name= "subs",
def=function(x,...)standardGeneric("subs")
)
#update: indel sign at polymorphic sites are reserved, feb 2018
#fixing a bug, May 2016
#fixing a bug, October 2016 (colnames(seq)<-names(polylist))
#subs method for Dna objects

setMethod(f="subs", signature= "Dna", definition=function (x,fifth=FALSE)
{
    
    seq<-toupper(x@sequence)
   
    seq[is.na(match(seq,c("A","T","C","G","-")))]<-"?"
    uniqs<-apply(seq,2,unique)
    if(is.matrix(uniqs))
    {
        uniqs<-lapply(seq_len(ncol(uniqs)), function(i) uniqs[,i])
        names(uniqs)<-colnames(seq)
    }
    whichpoly<-which(unlist(lapply(uniqs,length))>1)
    polymat<-seq[,whichpoly,drop=FALSE]
    
    if(ncol(polymat)>0)
    {
        
        if(!fifth)
        {
        polylist<- vector("list",0)
        
        for(i in 1:ncol(polymat))
        {
            
            substit<-unique(polymat[,i,drop=FALSE][polymat[,i,drop=FALSE]!="?"&polymat[,i,drop=FALSE]!="-"])
            
            if(length(substit)>1)
            {
                
                poly<-list(substit)
                names(poly)<-whichpoly[i]
                polylist<-c(polylist,poly)
            }
            
        }
        seq<-seq[,as.numeric(names(polylist)),drop=FALSE]
        colnames(seq)<-names(polylist)
        return(list(subsmat=seq,subs=polylist,subsmnum=length(polylist)))
        }
        
        if(fifth)
        {
            polylist<- vector("list",0)
            
            for(i in 1:ncol(polymat))
            {
                
                substit<-unique(polymat[,i,drop=FALSE][polymat[,i,drop=FALSE]!="?"])
                
                if(length(substit)>1)
                {
                    
                    poly<-list(substit)
                    names(poly)<-whichpoly[i]
                    polylist<-c(polylist,poly)
                }
                
            }
            seq<-seq[,as.numeric(names(polylist)),drop=FALSE]
            colnames(seq)<-names(polylist)
            return(list(subsmat=seq,subs=polylist,subsmnum=length(polylist)))
        }
        
        
        
    } else {
        return(list(subsmat=polymat,subs=list(),subsmnum=0))
    }
}
)





# internal function: fillendgaps

fillendgaps<-function(x,find="-",replace="?")
{
	nc<-ncol(x)
	for(i in 1:nrow(x))
	{
		f<-which(x[i,]==find)
		fl<-length(f)
		if(fl>0) 
		{
			if(fl==1)
			{
				beg<-f
				end<-f
				indell<-1
			} else 
			{
				fpoz<-c(2,f[2:fl]-f[1:(fl-1)])
				beg<-f [which(fpoz!=1)]
				end<-c(f[which(fpoz[-1]!=1)],f[fl])
			}
			
			fl<-length(beg)
			
			if(end[fl]==nc) x[i,beg[fl]:end[fl]]<- replace
			
		}
	}
	
	x
	
}

#internal function: alltest

alltest<-function(x,char="-") all(x==char)


#bugs fixed, feb 2018
# internal function: simple indel coding

sic<-function(x)
{
	
	allindel<-apply(x,2,alltest,char=5)
	x[,allindel]<-0
	
	indelmatrix<-matrix(NA,0,3)
	
	for(i in 1:nrow(x))
	{
        
		indels<-which(x[i,]==5)
		if(length(indels)>0) 
		{
			if(length(indels)==1)
			{
				beg<-indels
				end<-indels
				indell<-1
			} else 
			{
				indelpoz<-c(2,indels[2:length(indels)]-indels[1:(length(indels)-1)])
				beg<-indels [which(indelpoz!=1)]
				end<-c(indels[which(indelpoz[-1]!=1)],indels[length(indels)])
				indell<-1+end-beg
			}
            
            missdata<-x[i,]==0
            if(any(missdata))
            {
                wm<-which(missdata)
                st<-beg-1
                st[st<1]<-1
                en<-end+1
                en[en>ncol(x)]<-ncol(x)
                cakisma<-vector("logical",0)
                for(j in 1:length(beg))
                {
                    
                   
                    interv<-st[j]:en[j]
                    ck<-any(wm%in%interv)
                
                    if(ck) x[i,interv]<-0
                    cakisma<-c(cakisma,ck)
                    
                }
                beg<-beg[!cakisma]
                end<-end[!cakisma]
                indell<-indell[!cakisma]
            }
			
			indelmatrix<-rbind(indelmatrix,cbind(beg,end,indell))
			indelmatrix<-unique(indelmatrix)
			
		}
	}
	
	if(nrow(indelmatrix)==0) return(list(indels=indelmatrix,codematrix=matrix(NA,nrow(x),0)))

	indelmatrix<-indelmatrix[order(indelmatrix[,1],indelmatrix[,2]),,drop=FALSE]
	rownames(indelmatrix)<-1:nrow(indelmatrix)
	
	codematrix<-matrix(0,nrow(x),nrow(indelmatrix))
	element<-vector("list",nrow(indelmatrix))
	
	for(i in 1:ncol(codematrix))
	{
		beg<-indelmatrix[i,1]
		end<-indelmatrix[i,2]
		begother<-indelmatrix[-i,1,drop=TRUE]
		endother<-indelmatrix[-i,2,drop=TRUE]
        wi<-which(beg>=begother&end<=endother)
        if(length(wi)>0) if(is.null(names(wi))) names(wi)<-1
		element[[i]]<-as.numeric(names(wi))
		
	}	
	
	
	for(i in 1:ncol(codematrix))
	{
		beg<-indelmatrix[i,1]
		end<-indelmatrix[i,2]
		indelpart<-x[,beg:end]
		if(indelmatrix[i,3]==1) indelpart<-matrix(indelpart,,1)
		indels<- as.numeric(apply(indelpart,1,alltest,char=5))
        indels[apply(indelpart,1,alltest,char=0)]<--1
        codematrix[,i]<-indels
		
		el<-element[[i]]
		if(length(el)>0)
		{
			for(j in el)
			{
				beg<-indelmatrix[j,1]
				end<-indelmatrix[j,2]
				indelpart<-x[,beg:end]
				if(indelmatrix[j,3]==1) indelpart<-matrix(indelpart,,1)
				longerindels<- apply(indelpart,1,alltest,char=5)
				codematrix[longerindels,i]<--1
			}
			
		}
		
	}
	
	return(list(indels=indelmatrix,codematrix=codematrix))
	
}



#Generic indelcoder (only simple indel coding is available)

setGeneric (
name= "indelcoder",
def=function(x,...)standardGeneric("indelcoder")
)


#update: rownames given to codematrix, feb 2018.
#indelcoder method for Dna objects

setMethod(f="indelcoder", signature= "Dna", definition=function(x)
{
    
    x<-fillendgaps(x)
    seq<-as.numeric(x)
    sc<-sic(seq)
    rownames(sc$codematrix)<-names(x)
    return(list(indels=sc$indels,codematrix=sc$codematrix))
}
)



#Generic distance


setGeneric (
name= "distance",
def=function(x,...)standardGeneric("distance")
)

#updated for new subs method, feb 2018
#distance method for Dna objects

setMethod(f="distance", signature= "Dna", definition=function(x,subset=NULL,indels="sic")
{
	
	if(any(subset<1)|any(subset>nrow(x))) stop(paste("elements of 'subset' must be integers in the range [1,",nrow(x),"]",sep=""))
	
	indmet<- c("sic","5th","missing")
	matchmet <- pmatch(indels, indmet)
	if (is.na(matchmet)) 
	stop("invalid indel coding method", paste("", indels))
	indels<-indmet[matchmet]
	
	d<-nrow(x)
	if(nrow(x)==1) return(dist(0))
	dmat<-matrix(NA,d,d)
	rownames(dmat)<-x@seqnames
	
	
	comb<-combn(d,2)
	
	if(!is.null(subset))
	{
		newcomb<-comb
		newcomb[!is.na(match(newcomb,subset))]<-0
		comb<-comb[,newcomb[1,]*newcomb[2,]==0]
	}
	
	if(indels=="sic")
	{
		
		indeldat<-indelcoder(x)$codematrix
		
		if(ncol(indeldat)>0)
		{
			subsdat<-subs(x,fifth=FALSE)$subsmat
			
			for(i in 1:ncol(comb))
			{
				ind<-comb[,i]
				indeller<-indeldat[ind,,drop=FALSE]
				indeller<-indeller[,indeller[1,]!=-1&indeller[2,]!=-1,drop=FALSE]
				
				subslar<-subsdat[ind,,drop=FALSE]
                subslar<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!="-"&subslar[2,]!="-",drop=FALSE]
				
				dmat[ind[2],ind[1]]<-sum(indeller[1,]!=indeller[2,])+sum(subslar[1,]!=subslar[2,])
				
			}
		} else 
		{
			indels<-"missing"
			
		}
	}
	
	
	if(indels=="5th")
	{
		x<-fillendgaps(x)
		subsdat<-subs(x,fifth=TRUE)$subsmat
		for(i in 1:ncol(comb))
		{
			ind<-comb[,i]
			
			subslar<-subsdat[ind,,drop=FALSE]
			subslar<-subslar[,subslar[1,]!="?"&subslar[2,]!="?",drop=FALSE]
			
			dmat[ind[2],ind[1]]<-sum(subslar[1,]!=subslar[2,])
			
		}
		
		
	}
	
	if(indels=="missing")
	{
		subsdat<-subs(x)$subsmat
		for(i in 1:ncol(comb))
		{
			ind<-comb[,i]
			
			subslar<-subsdat[ind,,drop=FALSE]
			subslar<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!="-"&subslar[2,]!="-",drop=FALSE]
            
			dmat[ind[2],ind[1]]<-sum(subslar[1,]!=subslar[2,])
			
		}
		
	}
	
	
	return(as.dist(dmat))
}
)



#Generic polymorp

setGeneric (
name= "polymorp",
def=function(x,...)standardGeneric("polymorp")
)


#updated for new subs method, feb 2018
#polymorp method for Dna objects


setMethod(f="polymorp", signature= "Dna", definition=function(x,pair,indels="sic")
{
	if(nrow(x)<2) stop("at least two DNA sequences are required") 
	if(length(pair)!=2) stop("'pair' must be vector of length 2")
	if(any(pair<1)|any(pair>nrow(x))) stop(paste("elements of 'pair' must be integers in the range [1,",nrow(x),"]",sep=""))
	
	indmet<- c("sic","5th","missing")
	matchmet <- pmatch(indels, indmet)
	if (is.na(matchmet)) 
	stop("invalid indel coding method", paste("", indels))
	indels<-indmet[matchmet]
	
	seqmat<-x@sequence
	d<-nrow(seqmat)
	
	polylist<-list(list(),list())
	
	names(polylist)<-c("indels","subst")
	
	if(indels=="sic")
	{
		
		indel.obj<-indelcoder(x)
		indeldat<-indel.obj$codematrix
		subs.obj<-subs(x)
		subsdat<-subs.obj$subsmat
		
		indeller<-indeldat[pair,,drop=FALSE]
		
		nonmis<-which(indeller[1,]!=-1&indeller[2,]!=-1)
		indeller<-indeller[,nonmis,drop=FALSE]
		diff<-which(indeller[1,]!=indeller[2,])
		polyind<-nonmis[diff]
		ins<-indel.obj$indels[polyind,,drop=FALSE]
		indellist<-vector("list",0)
		if(nrow(ins)>0)
		{
			indellist<-vector("list",nrow(ins))
			names(indellist)<-ins[,1]
			for(i in 1:nrow(ins)) indellist[[i]]<-  seqmat[pair,ins[i,1]:ins[i,2],drop=FALSE]
			polylist[[1]]<-indellist
		}
		
		subslar<-subsdat[pair,,drop=FALSE]
		polysubs<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!="-"&subslar[2,]!="-"&subslar[1,]!=subslar[2,],drop=FALSE]
	
		subsmat<-polysubs
		subslist<-vector("list",0)
		if(ncol(subsmat)>0)
		{
			subslist<-vector("list",ncol(subsmat))
			names(subslist)<-colnames(subsmat)
			for(i in 1:ncol(subsmat)) subslist[[i]]<- subsmat[,i,drop=FALSE]	
		}
		polylist[[2]]<-subslist
		
		
	}
	
	
	if(indels=="5th")
	{
		
		x<-fillendgaps(x)
		subs.obj<-subs(x,fifth=TRUE)
		subsdat<-subs.obj$subsmat
		subslar<-subsdat[pair,,drop=FALSE]
		polysubs<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!=subslar[2,],drop=FALSE]
		subsmat<-polysubs		
		subslist<-vector("list",0)
		if(ncol(subsmat)>0)
		{
			
			subslist<-vector("list",ncol(subsmat))
			names(subslist)<-colnames(subsmat)
			for(i in 1:ncol(subsmat)) subslist[[i]]<- subsmat[,i,drop=FALSE]	
		}	
		polylist[[2]]<-subslist
	}
	
	if(indels=="missing")
	{
		subs.obj<-subs(x)
		subsdat<-subs.obj$subsmat
		subslar<-subsdat[pair,,drop=FALSE]
		subslar<-subsdat[pair,,drop=FALSE]
        polysubs<-subslar[,subslar[1,]!="?"&subslar[2,]!="?"&subslar[1,]!="-"&subslar[2,]!="-"&subslar[1,]!=subslar[2,],drop=FALSE]

		subsmat<-polysubs
	
		subslist<-vector("list",0)
		if(ncol(subsmat)>0)
		{
			
			subslist<-vector("list",ncol(subsmat))
			names(subslist)<-colnames(subsmat)
			for(i in 1:ncol(subsmat)) subslist[[i]]<- subsmat[,i,drop=FALSE]	
		}	
		polylist[[2]]<-subslist
	}
	
	
	
	return(polylist)
}
)



##fixing a bug in append, July 2016
#append method for Dna objects 

setMethod(f="append", signature= "Dna", definition=function(x,values)
{
	seq<- rbind(x@sequence,values@sequence)
	dnaobj<-new("Dna",sequence=seq,seqlengths=c(x@seqlengths,values@seqlengths),seqnames=c(x@seqnames,values@seqnames))
	dnaobj
})





#image method for Dna objects
setMethod(f="image", signature= "Dna", definition=function(x,all=FALSE,fifth=TRUE,col=c("#BFBFBF","#0B99FD","#FD0B0B","#11A808","#F5FD0B","#F8F8FF"),chars=TRUE,cex=1,show.names=TRUE,show.sites=TRUE,xlab="",ylab="",...)
{
    
   
    if(!all)
    {
        s<-subs(x,fifth=fifth)
        
        if(s$subsmnum==0)
        {
            
            warning("no polymorphism  found, all sites are displayed")
            all<-TRUE
            
        } else
        
        {
    
    s<-s$subsmat
    s<-s[nrow(s):1,,drop=FALSE]
    z<-as.numeric(as.dna(s))
    ele<-range(z)+1
    col<-col[ele[1]:ele[2]]
    
    if(show.names&&show.sites)
    {
        image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
        if(chars) text(col(z), row(z),label = s, cex=cex)
        mtext(rownames(s), side = 2, line = 0.1, at = 1:nrow(s),cex = cex, adj = 1, las = 1)
        mtext(colnames(s), side = 1, line = 0.1, at = 1:ncol(s),cex = cex, adj = 0.5, las = 1)
    } 
     if(show.names&&!show.sites)
    {
        image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
        if(chars) text(col(z), row(z),label = s, cex=cex)
       	mtext(rownames(s), side = 2, line = 0.1, at = 1:nrow(s),cex = cex, adj = 1, las = 1)
    
    
    }
     if(!show.names&&show.sites)
    {
        image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
        if(chars) text(col(z), row(z),label = s, cex=cex)
        mtext(colnames(s), side = 1, line = 0.1, at = 1:ncol(s),cex = cex, adj = 0.5, las = 1)
    
    
    }
      if(!show.names&&!show.sites)
    {
        image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
        if(chars) text(col(z), row(z),label = s, cex=cex)
    
    
    }
    
    
    
    
    }
    }
    
    
    if(all)
    {
        
        s<-as.matrix(x)
        s<-s[nrow(s):1,,drop=FALSE]
        z<-as.numeric(as.dna(s))
        ele<-range(z)+1
        col<-col[ele[1]:ele[2]]
         if(show.names&&show.sites)
        {
            image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
            if(chars) text(col(z), row(z),label = s, cex=cex)
            mtext(rownames(s), side = 2, line = 0.1, at = 1:nrow(x),cex = cex, adj =  1, las = 1)
            mtext(colnames(s), side = 1, line = 0.1, at = 1:ncol(x),cex = cex, adj = 0.5, las = 1)

        }       
        if(show.names&&!show.sites)
        {
            image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
            if(chars) text(col(z), row(z),label = s, cex=cex)
            mtext(rownames(s), side = 2, line = 0.1, at = 1:nrow(x),cex = cex, adj =  1, las = 1)

        }
        
               if(!show.names&&show.sites)
        {
            image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
            if(chars) text(col(z), row(z),label = s, cex=cex)
            mtext(colnames(s), side = 1, line = 0.1, at = 1:ncol(x),cex = cex, adj = 0.5, las = 1)

        }
        
               if(!show.names&&!show.sites)
        {
            image(x=1:ncol(z),y=1:nrow(z),z=t(z),col=col,yaxt="n",xaxt="n",xlab=xlab,ylab=ylab,...)
            if(chars) text(col(z), row(z),label = s, cex=cex)
       
        }
        
    }
    
    
    
}
)





#remove.gaps 
setGeneric (
name= "remove.gaps",
def=function(x,...)standardGeneric("remove.gaps")
)
#remove.gaps method for Dna objects

setMethod("remove.gaps","Dna", function(x,entire.col=FALSE)
{
    x<-toupper(x)
   
   if(entire.col)
   {
   	xl<-as.list(x)
    gp<-c()
    for(i in 1:nrow(x))
    {
        r<-xl[[i]]
        gaps<-which(r=="-")
        gp<-c(gaps,gp)
    }
        gp<-unique(gp)
    
    x[,-gp,as.matrix=FALSE]
    
   } else { 
   
   
    for(i in 1:nrow(x))
    {
        
        xi<-x[i,,as.matrix=TRUE]
        gaps<-which(xi=="-")
        x[i,]<- c(xi[-gaps,drop=FALSE],rep("-",length(gaps)))
        x@seqlengths[i]<-length(xi)-length(gaps)
    }

    
   x
   }
    
})







#Generic basecomp
setGeneric (
name= "basecomp",
def=function(x,...)standardGeneric("basecomp")
)
#base compositions of Dna objects
setMethod(f="basecomp", signature= "Dna", definition=function(x)
{
    seqnam<-names(x)
    x<-as.numeric(x)
    
    #total<- tabulate(as.numeric(x))[1:4]
    nseq<-nrow(x)
    M<-matrix(NA,nseq+1,4)
    
    for(i in 1:nseq)
    {
        tab<- tabulate(x[i,,drop=FALSE])[1:4]
        M[i,]<-tab/sum(tab)
        
    }
    
    tot<-tabulate(x)[1:4]
    M[(nseq+1),]<-tot/sum(tot)
    colnames(M)<-c("A","T","G","C")
    rownames(M)<-c(seqnam,"AVERAGE")
    M
    
}
)


#Generic function boot.dna
setGeneric (
name= "boot.dna",
def=function(x,...)standardGeneric("boot.dna")
)


#boot.dna method for Dna objects

setMethod("boot.dna","Dna", function(x,replacement=TRUE)
{

	l<-length(x)
	s<-sample(l,replace=TRUE)
    if(!replacement) s<-unique(s)
	s<-s[order(s)]
    as.dna(x@sequence[,s])
})



#new pairnei permutation support
#subset: n*2 matrix gives population pairs.

# internal function: pairwise nei raw D (D) : Nei s average number of differences between populations  (Nei and Li, 1979)


pairneidist<-function(x,populations,nperm=0,subset=NULL, showprogbar=FALSE)
{
    

    x<-as.matrix(x)
    pops<-unique(populations)
    npop<-length(pops)
    dmnam<-list(pops,pops)
    withinnei<-rep(NA,npop)
    if(is.null(subset)) sset<-1:npop else sset<-subset
	
  for(i in sset)
    {
        p<-populations==pops[i]
        xw<-x[p,p]
        n<-ncol(xw)
        if(!is.null(ncol(xw)))  nei <-mean(as.dist(xw)) else nei <- 0
        withinnei[i]<-nei
    }
 	
 	
 	
 
  
    neidistmat<-matrix(NA,npop,npop)
    
    comb<-combn(npop,2)
    if(!is.null(subset))
	{
		newcomb<-comb
		newcomb[!is.na(match(newcomb,subset))]<-0
		comb<-comb[,newcomb[1,]*newcomb[2,]==0]
	}
    
    if(!nperm)
    {
        
        pval<-NULL
        
        for(i in 1:ncol(comb))
        {
            
            p1<-populations==pops[comb[1,i]]
            p2<-populations==pops[comb[2,i]]
            nd<-mean(x=as.matrix(x[p1,p2]))
            neidistmat[comb[2,i],comb[1,i]]<-nd
            
            if(showprogbar) setTxtProgressBar(txtProgressBar(1,ncol(comb),style=3),i)
        }
        
        
    } else
    {
    
        pmat<-matrix(NA,npop,npop)
         
        for(i in 1:ncol(comb))
        {
                
               
                p1<-populations==pops[comb[1,i]]
                p2<-populations==pops[comb[2,i]]
                nd<-mean(x=as.matrix(x[p1,p2]))
                neidistmat[comb[2,i],comb[1,i]]<-nd
                
                p1p<-which(p1)
                p2p<-which(p2)
                pp<-c(p1p,p2p)
                l1<-length(p1p)
                
                pmat[comb[2,i],comb[1,i]]<-0
                for(j in 1:nperm)
                {
                    s<-sample(pp,replace=FALSE)
                    s1<-s[1:l1]
                    s2<-s[-c(1:l1)]
                    ndperm<-mean(x=as.matrix(x[s1,s2]))
                    pmat[comb[2,i],comb[1,i]]<- pmat[comb[2,i],comb[1,i]]+ as.numeric(ndperm>=nd)
                    
                }
                
                if(showprogbar) setTxtProgressBar(txtProgressBar(1,ncol(comb),style=3),i)
                
            }
           	
           
           	pval<-(pmat+1)/(nperm + 1)
            
            dimnames(pval)<-dmnam
           
                  }
    
 
    dimnames(neidistmat)<-dmnam
    
  	diag(neidistmat)<-withinnei
    listele<-list(neidist=neidistmat, p=pval)
    
    return(listele)
    
}


  
  
  
  
  
  
  
  
  
  
  
  



#new pairnei permutation support

#Generic pairnei
setGeneric (
name= "pairnei",
def=function(x,...)standardGeneric("pairnei")
)



#pairnei method for Dna object

setMethod("pairnei","Dna", function(x,populations,indels="sic",nperm=99,subset=NULL, showprogbar=FALSE)
{
    
    d<-distance(x,indels=indels)
    if(length(d)==0) stop("at least two DNA sequences are required")
    if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")
    if(length(populations)>nrow(x)) stop(paste("'populations' must be vector of length",nrow(x),"(equal or less than the number of DNA sequences)",sep=" "))
        if(any(subset<1)|any(subset>nrow(x))) stop(paste("elements of 'subset' must be integers in the range [1,",nrow(x),"]",sep=""))

    n<-pairneidist(d,populations,nperm=nperm, subset=subset, showprogbar= showprogbar)
    n
})


#pairnei method for matrix object

setMethod("pairnei","matrix", function(x,populations,nperm=99,subset=NULL,showprogbar=FALSE)
{
    if(nrow(x)==1) stop("dimensions of the distance matrix (x) must be greater than one")
    if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")
    if(length(populations)>nrow(x)) stop(paste("'populations' must be vector of length",nrow(x),"(equal or less than the dimensions of the distance matrix 'x')",sep=" "))
    if(any(subset<1)|any(subset>nrow(x))) stop(paste("elements of 'subset' must be integers in the range [1,",nrow(x),"]",sep=""))

    n<-pairneidist(as.dist(x),populations=populations,nperm=nperm, subset=subset,showprogbar= showprogbar)
    n
})


#pairnei method for dist object
#setOldClass("dist")

setMethod("pairnei","dist", function(x,populations,nperm=99,subset=NULL,showprogbar=FALSE)
{
    if(length(x)==0) stop("length of the distance object (x) must be greater than zero")
    if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")
    if(length(populations)>nrow(as.matrix(x))) stop(paste("'populations' must be vector of length",nrow(as.matrix(x)),"(equal or less than the dimensions of the distance matrix 'as.matrix(x)')",sep=" "))
        if(any(subset<1)|any(subset>nrow(as.matrix(x)))) stop(paste("elements of 'subset' must be integers in the range [1,",nrow(as.matrix(x)),"]",sep=""))

    n<-pairneidist(x,populations=populations,nperm=nperm, subset=subset,showprogbar= showprogbar)
    n	
})





#comb: n*2 matrix gives population pairs.
#pairfst
#using pegas amova and ade4 phiST internal


pairPhiSTdist<-function(x,populations,nperm=0,negatives=FALSE,subset=NULL,showprogbar=TRUE)
{
    
    ## Statphi is taken from package ade4 version 1.7-8 without any modification, author Sandrine Pavoine

    Statphi<-function(sigma) {
        f <- rep(0, length(sigma) - 1)
        if (length(sigma) == 3) {
            f <- rep(0, 1)
        }
        f[1] <- (sigma[length(sigma)] - sigma[length(sigma) -
        1])/sigma[length(sigma)]
        if (length(f) > 1) {
            s1 <- cumsum(sigma[(length(sigma) - 1):2])[-1]
            s2 <- sigma[(length(sigma) - 2):2]
            f[length(f)] <- sigma[1]/sigma[length(sigma)]
            f[2:(length(f) - 1)] <- s2/s1
        }
        return(f)
    }
    
    
    x<-as.matrix(x)
    
    pops<-unique(populations)
    npop<-length(pops)
    
    
    phiSTdistmat<-matrix(NA,npop,npop)
    
    
    
    
    pmat<-matrix(NA,npop,npop)
    
    dmnam<-list(pops,pops)
    dimnames(pmat)<-dmnam
    dimnames(phiSTdistmat)<-dmnam
    
    
    
    comb<-combn(npop,2)
    if(!is.null(subset))
	{
		newcomb<-comb
		newcomb[!is.na(match(newcomb,subset))]<-0
		comb<-comb[,newcomb[1,]*newcomb[2,]==0,drop=FALSE]
	}


    
    for(i in 1:ncol(comb))
    {
        
        
        p1<-which(populations==pops[comb[1,i]])
        p2<-which(populations==pops[comb[2,i]])
        
        pp<-c(p1,p2)
        
        distobj<-as.dist(x[pp,pp])
        dpopu<-as.factor(as.character(populations)[pp])
        
        
        if(length(dpopu)>2)
        {
            
            #use this temporary function until next update of the package pegas
            amo<-pegas.amova(distobj~dpopu,nperm=nperm,is.squared=TRUE)
            
             # amo<-pegas::amova(distobj~dpopu,nperm=nperm,is.squared=TRUE)
            if(nperm>0)
            {
                pmat[comb[2,i],comb[1,i]]<-amo$varcomp[1,2]
                sv<-sum(amo$varcomp[,1])
                if(sv>0) PhiST<-Statphi(c(amo$varcomp[,1],sv)) else PhiST=0
            } else {
                
                sv<-sum(amo$varcomp)
                if(sv>0) PhiST<-Statphi(c(amo$varcomp,sv)) else PhiST=0
            }
            
            phiSTdistmat[comb[2,i],comb[1,i]]<-PhiST
            
            
        } else
        {
            
            if(distobj>0) phiSTdistmat[comb[2,i],comb[1,i]]<-1
            if(distobj==0) phiSTdistmat[comb[2,i],comb[1,i]]<-0
            if(nperm>0) pmat[comb[2,i],comb[1,i]]<-1
            
        }
        
        if(showprogbar) setTxtProgressBar(txtProgressBar(1,ncol(comb),style=3),i)
        
    }

    if(nperm==0) pmat<-NULL
    
    if(!negatives) phiSTdistmat[phiSTdistmat<0]<-0
    listele<-list(PhiST=phiSTdistmat,p=pmat)
    
    return(listele)
    
}






#pairPhiST

#Generic pairPhiST
setGeneric (
name= "pairPhiST",
def=function(x,...)standardGeneric("pairPhiST")
)


#pairPhiST method for Dna object

setMethod("pairPhiST","Dna", function(x,populations,indels="sic",nperm=99,negatives=FALSE,subset=NULL,showprogbar=TRUE)
{
    
    d<-distance(x,indels=indels)
    if(length(d)==0) stop("at least two DNA sequences are required")
    if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")
    if(length(populations)>nrow(x)) stop(paste("'populations' must be vector of length",nrow(x),"(equal or less than the number of DNA sequences)",sep=" "))
    
    n<-pairPhiSTdist(d,populations,nperm=nperm,negatives=negatives,subset=subset,showprogbar=showprogbar)
        n
    
}
)



#pairPhiST method for matrix object

setMethod("pairPhiST","matrix", function(x,populations,nperm=99,negatives=FALSE,subset=NULL,showprogbar=TRUE)
{
    
    if(nrow(x)==1) stop("dimensions of the distance matrix (x) must be greater than one")
    if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")
    if(length(populations)>nrow(x)) stop(paste("'populations' must be vector of length",nrow(x),"(equal or less than the dimensions of the distance matrix 'x')",sep=" "))
    
    n<-pairPhiSTdist(as.dist(x),populations=populations,nperm=nperm,negatives=negatives, subset = subset,showprogbar=showprogbar)
    n
})


#pairPhiST method for dist object

setMethod("pairPhiST","dist", function(x,populations,nperm=99,negatives=FALSE, subset =NULL,showprogbar=TRUE)
{
    if(length(x)==0) stop("length of the distance object (x) must be greater than zero")
    if(length(populations)==1) stop("'populations' must be vector of length >1 (at least two populations are required)")
    if(length(populations)>nrow(as.matrix(x))) stop(paste("'populations' must be vector of length",nrow(as.matrix(x)),"(equal or less than the dimensions of the distance matrix 'as.matrix(x)')",sep=" "))
   
    n<-pairPhiSTdist(x,populations=populations,nperm=nperm,negatives=negatives, subset = subset,showprogbar=showprogbar)
    n
})

#new function
#internal function indel method test
indtest<-function(indels="sic")
{
indmet<- c("sic","5th","missing")
matchmet <- pmatch(indels, indmet)
if (is.na(matchmet)) 
stop("invalid indel coding method", paste("", indels))
indels<-indmet[matchmet]
indels
}




#Coerce Dna object to a DNAbin object
setMethod(f="as.DNAbin", signature= "Dna", definition=function(x,endgaps=TRUE)
{
	x<-tolower(x)
	rng<-range(x)
	if(rng[1]==rng[2]|endgaps)
	{
		x<-x@sequence
	} else {
		x<-as.list(x)
	}
	ape::as.DNAbin(x)
}
)


#new function
#internal function contr
contr<-function(indels="sic")
{
    indels<-indtest(indels)
    
    if(indels=="sic")
    {
        
        contrast<-matrix(data = c(1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1,
        1,1,1,1,1,1,
        1,1,1,1,1,1,
        1,1,1,1,1,1),
        ncol = 6, byrow = TRUE)
        
        dimnames(contrast)<-list(c("A","C","G","T","0","1","-","?","-1"), c("A","C","G","T","0","1"))
    }
    if(indels=="5th")
    {
        
        contrast<-matrix(data = c(1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1,
        1,1,1,1,1),
        ncol = 5, byrow = TRUE)
        
        dimnames(contrast)<-list(c("A","C","G","T","-","?"), c("A","C","G","T","-"))
        
        
    }
    
    
    if(indels=="missing")
    {
        
        contrast<-matrix(data = c(1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1,
        1,1,1,1,
        1,1,1,1),
        ncol = 4, byrow = TRUE)
        
        dimnames(contrast)<-list(c("A","C","G","T","-","?"), c("A","C","G","T"))
        
    }
    
    contrast
}




#as.phydat method for Dna objects
#setOldClass("as.phydat")

setMethod(f="as.phyDat", signature= c("Dna"), definition=function(x,indels="sic",...)
{

x<-toupper(x)
indels<-indtest(indels)

seq<-x@sequence

  
if(indels=="sic")
  {
		inds<-indelcoder(x)$codematrix
            contrast<-contr(indels)
            seqind<-cbind(seq, inds)
            xphy<-phangorn::as.phyDat(seqind, type="USER",contrast= contrast,...)
                
    } 
if(indels=="5th")
    {
		contrast<-contr(indels)
            xphy<-phangorn::as.phyDat(seq, type="USER",contrast= contrast,...)
    }
if(indels=="missing")
    {
		contrast<-contr(indels)
            xphy<-phangorn::as.phyDat(seq, type="DNA",...)
    }


    return(xphy)
    
    

}
)





#update amova
#fix pegas::amova, use this temporary function until next update of pegas
#Feb 28, 2020

#temporary internal function pegas.amova
pegas.amova<-function (formula, data = NULL, nperm = 1000, is.squared = FALSE)
{
    y.nms <- as.character(as.expression(formula[[2]]))
    rhs <- formula[[3]]
    gr.nms <- as.character(as.expression(rhs))
    if (length(rhs) > 1)
    gr.nms <- unlist(strsplit(gr.nms, "/"))
    data.env <- if (is.null(data)) environment(formula) else as.environment(data)
    if (any(sapply(gr.nms, function(x) !is.factor(get(x, envir = data.env)))))
    warning("elements in the rhs of the formula are not all factors")
    gr <- as.data.frame(sapply(gr.nms, get, envir = data.env),stringsAsFactors = TRUE)
    y <- get(y.nms, envir = environment(formula))
    if (any(is.na(y)))
    warning("at least one missing value in the distance object.")
    if (!is.squared)
    y <- y^2
    if (class(y) == "dist")
    y <- as.matrix(y)
    if (!is.matrix(y))
    stop("the lhs of the formula must be either a matrix or an object of class 'dist'.")
    n <- dim(y)[1]
    Nlv <- length(gr)
    if (nperm) {
        o <- do.call("order", gr)
        gr <- gr[o, , drop = FALSE]
        if (Nlv > 1) {
            f <- function(x) {
                lp <- as.character(unique(x))
                factor(match(as.character(x), lp))
            }
            for (i in 2:Nlv) gr[, i] <- f(gr[, i])
        }
        y <- y[o, o]
    }
    foo <- function(x, index) unlist(lapply(split(x, index),
    sum))
    getSSD <- function(y, gr, Nlv, N, n) {
        SSD <- numeric(Nlv + 2)
        SSD[Nlv + 2] <- sum(y/(2 * n))
        for (i in 1:Nlv) {
            p <- gr[, i]
            SSD[i + 1] <- sum((y/(2 * N[[i]])[p])[outer(p, p,
            "==")])
        }
        if (Nlv > 1)
        for (i in 2:Nlv) SSD[i] <- SSD[i] - SSD[i + 1]
        SSD[1] <- SSD[Nlv + 2] - sum(SSD[-(Nlv + 2)])
        SSD
    }
    getDF <- function(gr, Nlv, N, n) {
        df <- numeric(Nlv + 2)
        df[1:Nlv] <- unlist(lapply(N, length))
        df[Nlv + 1] <- n
        for (i in (Nlv + 1):2) df[i] <- df[i] - df[i - 1]
        df[1] <- df[1] - 1
        df[Nlv + 2] <- n - 1
        df
    }
    getNcoefficient <- function(gr, Nlv, N, n) {
        Nig <- N[[Nlv]]
        Nig2 <- Nig^2
        npop <- length(Nig)
        if (Nlv == 1)
        ncoef <- (n - sum(Nig2)/n)/(npop - 1)
        else {
            if (Nlv == 2) {
                ncoef <- numeric(3)
                G <- nlevels(gr[, 1])
                g <- gr[, 1][match(1:npop, as.integer(gr[, 2]))]
                npopBYgr <- tabulate(g)
                A <- sum(foo(Nig2, g)/foo(Nig, g))
                ncoef[1] <- (n - A)/sum(npopBYgr - 1)
                ncoef[2] <- (A - sum(Nig2)/n)/(G - 1)
                ncoef[3] <- (n - sum(foo(Nig, g)^2/n))/(G - 1)
            }
            else {
                ncoef <- numeric(Nlv + 1)
                ncoef[Nlv] <- (n - sum(Nig2)/n)/(npop - 1)
                ncoef[Nlv + 1] <- 1
                for (i in 1:(Nlv - 1)) {
                    group <- gr[, i]
                    g <- group[match(1:npop, as.integer(gr[, i +
                    1]))]
                    A <- sum(foo(Nig, g)^2)/sum(foo(Nig, g))
                    ncoef[i] <- (n - A)/(nlevels(group) - 1)
                }
            }
        }
        names(ncoef) <- letters[1:length(ncoef)]
        ncoef
    }
    getVarComp <- function(MSD, Nlv, ncoef) {
        if (Nlv == 1)
        sigma2 <- c((MSD[1] - MSD[2])/ncoef, MSD[2])
        else {
            sigma2 <- numeric(Nlv + 1)
            if (Nlv == 2) {
                sigma2[3] <- MSD[3]
                sigma2[2] <- (MSD[2] - sigma2[3])/ncoef[1]
                sigma2[1] <- (MSD[1] - MSD[3] - ncoef[2] * sigma2[2])/ncoef[3]
            }
            else {
                sigma2[Nlv + 1] <- MSD[Nlv + 1]
                for (i in Nlv:1) {
                    sel <- i:(Nlv + 1)
                    sigma2[i] <- (MSD[i] - sum(ncoef[sel] * sigma2[sel]))/ncoef[i]
                }
            }
        }
        names(sigma2) <- c(names(gr), "Error")
        sigma2
    }
    N <- lapply(gr, tabulate)
    SSD <- getSSD(y, gr, Nlv, N, n)
    df <- getDF(gr, Nlv, N, n)
    MSD <- SSD/df
    ncoef <- getNcoefficient(gr, Nlv, N, n)
    sigma2 <- getVarComp(MSD, Nlv, ncoef)
    res <- list(tab = data.frame(SSD = SSD, MSD = MSD, df = df,
    row.names = c(names(gr), "Error", "Total"),stringsAsFactors = FALSE), varcoef = ncoef,
    varcomp = sigma2, call = match.call())
    class(res) <- "amova"
    if (nperm) {
        rSigma2 <- matrix(0, nperm, length(sigma2))
        j <- if (Nlv == 1)
        1:2
        else Nlv + 1
        for (i in 1:nperm) {
            rY <- ape::perm.rowscols(y, n)
            rSSD <- getSSD(rY, gr, Nlv, N, n)
            rSigma2[i, j] <- getVarComp(rSSD/df, Nlv, ncoef)[j]
        }
        if (Nlv > 1) {
            j <- Nlv
            L <- lapply(levels(gr[, j - 1]), function(x) which(gr[,
            j - 1] == x))
            for (i in 1:nperm) {
                rind <- unlist(lapply(L, sample))
                rY <- y[rind, rind]
                rSSD <- getSSD(rY, gr, Nlv, N, n)
                rSigma2[i, j] <- getVarComp(rSSD/df, Nlv, ncoef)[j]
            }
            if (Nlv > 2) {
                for (j in (Nlv - 1):2) {
                    above <- gr[, j - 1]
                    L <- lapply(levels(above), function(x) which(above ==
                    x))
                    for (i in 1:nperm) {
                        rind <- integer(0)
                        for (k in L) rind <- c(rind, sample(k))
                        rind <- unlist(lapply(L, sample))
                        rY <- y[rind, rind]
                        rGR <- gr[rind, ]
                        rN <- lapply(rGR, tabulate)
                        rSSD <- getSSD(rY, rGR, Nlv, rN, n)
                        rDF <- getDF(rGR, Nlv, rN, n)
                        rNcoef <- getNcoefficient(rGR, Nlv, rN, n)
                        rSigma2[i, j] <- getVarComp(rSSD/rDF, Nlv,
                        rNcoef)[j]
                    }
                }
            }
            N2 <- N[[2]]
            Higr <- gr[, 1][cumsum(N2)]
            rGR <- gr
            for (i in 1:nperm) {
                rGR[, 1] <- unlist(mapply(rep, sample(Higr),
                each = N2, SIMPLIFY = FALSE))
                rN <- lapply(rGR, tabulate)
                rSSD <- getSSD(rY, rGR, Nlv, rN, n)
                rDF <- getDF(rGR, Nlv, rN, n)
                rNcoef <- getNcoefficient(rGR, Nlv, rN, n)
                rSigma2[i, 1] <- getVarComp(rSSD/rDF, Nlv, rNcoef)[1]
            }
        }
        P <- numeric(Nlv + 1)
        for (j in 1:(Nlv + 1)) P[j] <- sum(rSigma2[, j] >= sigma2[j])/(nperm +
        1)
        P[Nlv + 1] <- NA
        res$varcomp <- data.frame(sigma2 = res$varcomp, P.value = P,stringsAsFactors = FALSE)
    }
    res
}

































