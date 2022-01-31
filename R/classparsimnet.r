### class Parsimnet, functions and methods
##march 2015, cnrakt
##function .TempletonProb is taken from package pegas version 0.6, authors Emmanuel Paradis, Klaus Schliep.
##Some internal structures of pieplot is taken from package network with modifications, author Carter T. Butts.

#################May 2016, cnrakt#################

#fixing a bug/update, May 2016
#internal function: calculate steplimit

steplimit<-function(seqlength,prob=.95)
{
    
    # .TempletonProb is taken from package pegas version 0.6, authors Emmanuel Paradis, Klaus Schliep.
    .TempletonProb <- function(j, S, b = 2, r = 1)
    {
        br <- b * r
        P <- numeric(max(j))
        L_jm <- function(q, j, m) {
            jm1 <- j - 1
            qonbr <- q/br
            (2*q)^jm1 * (1 - q)^(2*m + 1) * (1 - qonbr) *
            (2 - q*(br + 1)/br)^jm1 *
            (1 - 2*q*(1 - qonbr))
        }
        for (i in seq_along(P)) {
            M <- S - i
            denom <- integrate(L_jm, 0, 1, j = i, m = M)$value
            ## eq.7 from Templeton et al. 1992:
            out <- integrate(function(q) q*L_jm(q, j = i, m = M), 0, 1)$value/denom
            P[i] <- 1 - out
        }
        cumprod(P)[j]
    }
    
    i<-1
    test<-1
    testvec<-NULL
    while(test>=prob&&i<seqlength)
    {
        test<-.TempletonProb(i,seqlength)
        if(is.nan(test)) break()
        i<-i+1
        testvec<-c(testvec,test)
        
    }
    testvec<-testvec[-(i-1)]
    return(testvec)
}


#internal function: tryconnect

tryconnect<-function(dmat,newhapnam,step,prevp,nexp,prevpair,nhap)
{
    
    if(step>1)
    {
        
        dmat<- fillna(dmat,prevpair,prevp,nexp,step)
        temproxhapsobj<-xhaps(dmat,prevpair,prevp,nexp,step,newhapnam,rownames=c(prevp,nexp))
        temprodmat<-temproxhapsobj$d
        temproprevpair<-temproxhapsobj$prevpair
        tempronewhapnam<-nrow(temprodmat)
        
        
        net<-list(d=temprodmat,newhapnam=tempronewhapnam,step=step,prevp=prevp,nexp=nexp,prevpair=temproprevpair)
        
    } else {
        
        net<-list(d=dmat,newhapnam=newhapnam,step=step,prevp=prevp,nexp=nexp,prevpair=prevpair)
        
        
    }
    
    
    return(net)
    
}


#internal function: fillna

fillna<-function(dmat,prevpair,prevp,nexp,step)
{
    nas<-is.na(dmat[nexp,])
    nas[nexp]<-FALSE
    filna<-dmat[prevp,nas]
    filna[is.na(filna)]<-0
    dmat[nas,nexp]<-dmat[nexp,nas]<-filna+step
    
    return(dmat)
    
}

#bug fixed in calculating , May 2016
#internal function: xhap

xhaps<-function(dmat,prevpair,prevp,nexp,step,newhapnam,rownames)

{
    
    prevpair<-c(prevpair,nexp)
    prevpair<-prevpair[order(prevpair)]
    
    convec<-dmat[prevp,]
    convec[nexp]<-step
    diffvec<-rep(1,length(convec))
    diffvec[prevp]<-1
    diffvec[nexp]<--1
    convec[prevp]<-0
    tempprevpair<-prevpair
    eskiprevpair<-prevpair
    rownam<-rownames
    
    
    komsi<-numeric(0)
    
    for(s in 1:(step-1))
    {
        
        if(s>1) komsi<-(s-1):1
        addp<-matrix(c(convec+diffvec*s,komsi),1,)
        navec<-matrix(rep(NA,nrow(dmat)),1,)
        addpm<-matrix(addp[,tempprevpair],1,)
        navec[,tempprevpair]<-addpm
        rownames(navec)<-paste(rownam[1],rownam[2],sep="_")
        dmattemp<-rbind(dmat,navec)
        diffmatrix<- t(-dmattemp[newhapnam+1,]+t(dmattemp[-c(newhapnam+1),]))
        reppoi<-apply(abs(diffmatrix),1,sum,na.rm=TRUE)
        newhapnam<-newhapnam+1
        rownames(addp)<-paste(rownam[1],rownam[2],sep="_")
        rownames(addpm)<-rownames(addp)
        dmat<-rbind(dmat,navec)
        dmat<-cbind(dmat,t(cbind(navec,NA)))
        tempprevpair<-c(tempprevpair,newhapnam)
        
        
    }
    
    list(d=abs(dmat),prevpair=tempprevpair,newhapnam=newhapnam)
    
}



#update/fixing a bug
#internal function:calcclust

calcclust<-function(dmat,step=1)
{
    diag(dmat)<-0
    pairs<-which(dmat<=step,arr.ind=TRUE)
    nr<-nrow(dmat)
    komsilist<-vector("list",nr)
    
    for(i in 1:nr) komsilist[[i]]<-pairs[pairs[,2]==i,1]
    
    
    iter<-TRUE
    
    while(iter)
    {
        
        
        prevkomsilist<-komsilist
        
        for(i in 1:nr)
        {
            k<-komsilist[[i]]
            unk<-c()
            for(t in k) unk<-c(unk, komsilist[[t]])
            komsilist[[i]]<-unique(unk)
        }
        
        nr<-length(komsilist)
        for(i in nr:1)
        {
            k<-komsilist[[i]]
            unk<-c()
            for(t in k) unk<-c(unk, komsilist[[t]])
            komsilist[[i]]<-unique(unk)
        }
        
        
        iter<- !identical(prevkomsilist,komsilist)
        
    }
    
    komsilist<-unique(komsilist)
    l<-length(komsilist)
    for(i in 1:l) komsilist[[i]]<-sort(komsilist[[i]])
    komsilist<-unique(komsilist)
    l<-length(komsilist)
    
    mat<-matrix(NA,nr,2)
    for(i in 1: l)
    {
        mat[komsilist[[i]],1]<-i
        mat[komsilist[[i]],2]<-komsilist[[i]]
        
    }
    
    return(mat)
}




#internal function: recdistG

recdistG<-function(d,dmat,clustmat)
{
    
    nhap<-nrow(d)
    tempdmat<-dmat
    tempdmat[1:nhap,1:nhap]<-d
    clusters<-unique(clustmat[,1])
    
    for(i in clusters)
    {
        ind<-clustmat[clustmat[,1]==i,2]
        tempdmat[ind,ind]<-dmat[ind,ind]
    }
    return(tempdmat)
    
}


#fixing a bug in tryconnectG,  Recycling array of length 1 in vector-array arithmetic is deprecated warning.
#fixing a bug,  feb 2018
#internal function: tryconnectG
tryconnectG<-function(dmat,prevset,nexpset,nhap)
{
    scorevec<-vector("numeric",0)
    parsimscorevec<-vector("numeric",0)
    nexpvec<-vector("numeric",0)
    prevpvec<-vector("numeric",0)
    xvec<-vector("numeric",0)
    jumplist<-vector("list",0)
    
    
    for(i in 1:length(prevset))
    {
        
        prevp<-prevset[i]
        
        for(j in 1:length(nexpset))
        {
            
            nexp<-nexpset[j]
            nexpvec<-c(nexpvec,nexp)
            prevpvec<-c(prevpvec,prevp)
            
            if(prevp<=nhap)
            {
                if(nexp<=nhap)
                {
                    
                    x<-dmat[prevp,nexp]
                    jump<-dmat
                    
                    
                    jump[nexp,prevset]<-jump[prevset,nexp]<-jump[prevp,prevset]+x
                    jump[nexp,prevp]<-jump[prevp,nexp]<-x
                    if(length(nexpset)>1)
                    {
                        remmat<-matrix(jump[nexp,prevset],length(nexpset[-j]),length(prevset),byrow=TRUE)+jump[nexpset[-j],nexp]
                        jump[nexpset[-j],prevset]<-remmat
                        jump[prevset,nexpset[-j]]<-t(remmat)
                    }
                    test<-jump[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]-dmat[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]
                    score<-sum(test)
                    parssc<-x+sum(jump==1,na.rm=TRUE)/2
                    
                    if(any(test<0,na.rm=TRUE))
                    {
                        jump[nexpset, prevset]<-jump[prevset,nexpset]<-jump[nexpset, prevset]-min(test)
                        score<-sum(test-min(test))
                        x<-x-min(test)
                        parssc<-x+sum(jump==1,na.rm=TRUE)/2
                        
                        
                    }
                    parsimscorevec<-c(parsimscorevec,parssc)
                    scorevec<-c(scorevec,score)
                    xvec<-c(xvec,x)
                    jumplist<-c(jumplist,list(jump[nexpset, prevset]))
                }
                
                if(nexp>nhap)
                {
                    
                    minkomsi<-which.min(dmat[nexp,nexpset[nexpset<=nhap]])
                    komsi<-nexpset[nexpset<=nhap][minkomsi]
                    x<- dmat[prevp,komsi]-dmat[nexp,komsi]
                    jump<-dmat
                    jump[nexp,prevset]<-jump[prevset,nexp]<-jump[prevp,prevset]+x
                    
                    
                    jump[nexp,prevp]<-jump[prevp,nexp]<-x
                    if(length(nexpset)>1)
                    {
                        remmat<-matrix(jump[nexp,prevset],length(nexpset[-j]),length(prevset),byrow=TRUE)+jump[nexpset[-j],nexp]
                        jump[nexpset[-j],prevset]<-remmat
                        jump[prevset,nexpset[-j]]<-t(remmat)
                    }
                    test<-jump[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]-dmat[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]
                    score<-sum(test)
                    parssc<-x+sum(jump==1,na.rm=TRUE)/2
                    
                    if(any(test<0,na.rm=TRUE))
                    {
                        jump[nexpset, prevset]<-jump[prevset,nexpset]<-jump[nexpset, prevset]-min(test)
                        score<-sum(test-min(test))
                        x<-x-min(test)
                        parssc<-x+sum(jump==1,na.rm=TRUE)/2
                        
                        
                    }
                    parsimscorevec<-c(parsimscorevec,parssc)
                    scorevec<-c(scorevec,score)
                    xvec<-c(xvec,x)
                    jumplist<-c(jumplist,list(jump[nexpset, prevset]))
                }
            }
            
            
            if(prevp>nhap)
            {
                if(nexp<=nhap)
                {
                    minkomsi<-which.min(dmat[prevp,prevset[prevset<=nhap]])
                    komsi<-prevset[prevset<=nhap][minkomsi]
                    x<- dmat[nexp,komsi]-dmat[prevp,komsi]
                    jump<-dmat
                    jump[nexp,prevset]<-jump[prevset,nexp]<-jump[prevp,prevset]+x
                    
                    
                    
                    jump[nexp,prevp]<-jump[prevp,nexp]<-x
                    if(length(nexpset)>1)
                    {
                        remmat<-matrix(jump[nexp,prevset],length(nexpset[-j]),length(prevset),byrow=TRUE)+jump[nexpset[-j],nexp]
                        jump[nexpset[-j],prevset]<-remmat
                        jump[prevset,nexpset[-j]]<-t(remmat)
                    }
                    test<-jump[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]-dmat[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]
                    score<-sum(test)
                    parssc<-x+sum(jump==1,na.rm=TRUE)/2
                    
                    if(any(test<0,na.rm=TRUE))
                    {
                        jump[nexpset, prevset]<-jump[prevset,nexpset]<-jump[nexpset, prevset]-min(test)
                        score<-sum(test-min(test))
                        x<-x-min(test)
                        parssc<-x+sum(jump==1,na.rm=TRUE)/2
                        
                    }
                    parsimscorevec<-c(parsimscorevec,parssc)
                    scorevec<-c(scorevec,score)
                    xvec<-c(xvec,x)
                    jumplist<-c(jumplist,list(jump[nexpset, prevset]))
                }
                
                
                
                
                if(nexp>nhap)
                {
                    minkomsiP<-which.min(dmat[prevp,prevset[prevset<=nhap]])
                    komsiP<-prevset[prevset<=nhap][minkomsiP]
                    minkomsiN<-which.min(dmat[nexp,nexpset[nexpset<=nhap]])
                    komsiN<-nexpset[nexpset<=nhap][minkomsiN]
                    x<- dmat[komsiN,komsiP] - dmat[prevp,komsiP]-dmat[nexp,komsiN]
                    jump<-dmat
                    jump[nexp,prevset]<-jump[prevset,nexp]<-jump[prevp,prevset]+x
                    
                    
                    
                    jump[nexp,prevp]<-jump[prevp,nexp]<-x
                    if(length(nexpset)>1)
                    {
                        remmat<-matrix(jump[nexp,prevset],length(nexpset[-j]),length(prevset),byrow=TRUE)+jump[nexpset[-j],nexp]
                        jump[nexpset[-j],prevset]<-remmat
                        jump[prevset,nexpset[-j]]<-t(remmat)
                    }
                    test<-jump[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]-dmat[nexpset[nexpset<=nhap], prevset[prevset<=nhap]]
                    score<-sum(test)
                    parssc<-x+sum(jump==1,na.rm=TRUE)/2
                    
                    if(any(test<0,na.rm=TRUE))
                    {
                        jump[nexpset, prevset]<-jump[prevset,nexpset]<-jump[nexpset, prevset]-min(test)
                        score<-sum(test-min(test))
                        x<-x-min(test)
                        parssc<-x+sum(jump==1,na.rm=TRUE)/2
                        
                    }
                    parsimscorevec<-c(parsimscorevec,parssc)
                    scorevec<-c(scorevec,score)
                    xvec<-c(xvec,x)
                    jumplist<-c(jumplist,list(jump[nexpset, prevset]))
                    
                }
                
            }
            
            
        }
        
        
    }
    
    minnetlen<-min(parsimscorevec)
    
    maxpar<-which(parsimscorevec==minnetlen)
    best<-maxpar[which.min(scorevec[maxpar])]
    
    return(list(scorevec=scorevec[best],parsimvec=parsimscorevec[best],prevpvec=prevpvec[best],nexpvec=nexpvec[best],xvec=xvec[best],jump=jumplist[best]))
    
}


#bug fixed in calculating maximum connection limit,it was removed from clustersG function May 2016
#internal function: clustersG

clustersG<-function(d,nhap)
{
    test<-geodist.corr(d)-d
    test[is.na(test)]<-Inf
    if(length(test)>0) if(min(test,na.rm=TRUE)<0) stop("inconsistency in distance matrix, missing data in DNA sequences may cause problem when using algorithmic method. Consider using heuristic method.")
    
    dmat<-d
    clustmat<- calcclust(dmat,1)
    clusters<-unique(clustmat[,1])
    if(length(clusters)==1) loop<-FALSE else loop<-TRUE
    lind<-0
    while(loop)
    {
        
        
        
        dmat<-recdistG(d,dmat,clustmat)
        tempdmat<-dmat
        newhapnam<-nrow(dmat)
        clusters<-unique(clustmat[,1])
        clustcombn<-combn(clusters,2)
        scorevec<-vector("numeric",0)
        parsimvec<-vector("numeric",0)
        prevpvec<-vector("numeric",0)
        nexpvec<-vector("numeric",0)
        xvec<-vector("numeric",0)
        nexpsetlist<-vector("list",0)
        prevsetlist<-vector("list",0)
        jumplist<-vector("list",0)
        
        for(k in 1:ncol(clustcombn))
        
        {
            ind<-clustcombn[,k]
            size1<-sum(clustmat[,1]==ind[1])
            size2<-sum(clustmat[,1]==ind[2])
            if(size1>size2)
            {
                clA<-clustmat[clustmat[,1]==ind[2],2]
                clB<-clustmat[clustmat[,1]==ind[1],2]
            }
            
            if(size1<=size2)
            {
                clA<-clustmat[clustmat[,1]==ind[1],2]
                clB<-clustmat[clustmat[,1]==ind[2],2]
            }
            
            nexpset<-clA
            prevset<-clB
            
            trG.obj<-tryconnectG(dmat,prevset,nexpset,nhap)
            scorevec<-c(scorevec,trG.obj$scorevec)
            parsimvec<-c(parsimvec,trG.obj$parsimvec)
            prevpvec<-c(prevpvec,trG.obj$prevpvec)
            nexpvec<-c(nexpvec,trG.obj$nexpvec)
            xvec<-c(xvec,trG.obj$xvec)
            nexpsetlist<-c(nexpsetlist,list(nexpset))
            prevsetlist<-c(prevsetlist,list(prevset))
            jumplist<-c(jumplist,trG.obj$jump)
            
        }
        
        
        
        minnetlen<-min(parsimvec)
        maxpar<-which(parsimvec==minnetlen)
        best<-maxpar[which.min(scorevec[maxpar])]
        prevset<-prevsetlist[[best]]
        nexpset<-nexpsetlist[[best]]
        prevp<-prevpvec[best]
        nexp<-nexpvec[best]
        x<-xvec[best]
        jump<-jumplist[[best]]
        
        
        tempdmat[nexpset,prevset]<-jump
        tempdmat[prevset,nexpset]<-t(jump)
        
        
        if(x>1)
        {
            
            tempnet<-tryconnect(tempdmat,newhapnam,step=x,prevp=prevp,nexp,c(nexp,prevp),nhap)
            tempdmat<-tempnet$d
            newhaps<-(nrow(dmat)+1):nrow(tempdmat)
            
            remprevpset<-prevset[!prevset==prevp]
            
            if(length(remprevpset)>0)
            {
                
                addh<-tempdmat[prevp,remprevpset]
                addxh<-tempdmat[prevp,newhaps]
                recdmat<-matrix(addxh,length(addh),length(addxh),byrow=TRUE)+addh
                
                tempdmat[remprevpset,newhaps]<-recdmat
                tempdmat[newhaps,remprevpset]<-t(recdmat)
            }
            remnexpset<-nexpset[!nexpset==nexp]
            
            if(length(remnexpset)>0)
            {
                
                addh<-tempdmat[nexp,remnexpset]
                addxh<-tempdmat[nexp,newhaps]
                recdmat<-matrix(addxh,length(addh),length(addxh),byrow=TRUE)+addh
                
                tempdmat[remnexpset,newhaps]<-recdmat
                tempdmat[newhaps,remnexpset]<-t(recdmat)
            }
        }
        
        
        dmat<-tempdmat
        clustmat<- calcclust(dmat,1)
        
        
        if(length(unique(clustmat[,1]))==1) loop<-FALSE
        
        
    }
    
    return(dmat)
}




#bug fixed in calculating slot d, May 2016
#internal function: parsnet
#prob can be NULL
parsnet<-function(d,seqlength,prob=.95)
{
    
    d<-as.matrix(d)
    diag(d)<-NA
    rwnd<-rownames(d)
    if(is.null(rwnd)) rownames(d)<-1:nrow(d)
    nhaps<-vector("numeric",0)
    nhap<-nrow(d)
    if(!is.null(prob))
    {
        tempprobs<-steplimit(seqlength,prob)
        conlim<-length(tempprobs)
        clusters<-calcclust(d,step=conlim)
        
        if(all(clusters[,1]==1))
        {
            
            dmat<-clustersG(d,nhap)
            dmat<-geodist.corr(dmat)
            nm<-paste("net",1,sep="")
            grouplist<-list(1:nrow(d))
            dmatlist<-list(dmat)
            names(grouplist)<-nm
            names(dmatlist)<-nm
            numnet<-length(dmatlist)
            for(i in 1:numnet) rownames(dmatlist[[i]])<- ren.dupli(rownames(dmatlist[[i]]))
            
            return(list(d=dmatlist,tempProbs=tempprobs,conlimit=conlim,prob=prob,nhap=nhap,rowindex=grouplist))
        } else {
            
            cl<-unique(clusters[,1])
            dmatlist<-vector("list",length(cl))
            grouplist<-vector("list",length(cl))
            for(i in 1: length(cl))
            {
                group<-clusters[clusters[,1]==cl[i],2]
                dma<-d[group,group,drop=FALSE]
                dmat<-clustersG(dma,length(group))
                dmat<-geodist.corr(dmat)
                grouplist[[i]]<-group
                dmatlist[[i]]<-dmat
                nhaps<-c(nhaps,length(group[group<=nhap]))
            }
            nm<-paste("net",1:length(cl),sep="")
            names(grouplist)<-nm
            names(dmatlist)<-nm
            numnet<-length(dmatlist)
            for(i in 1:numnet) rownames(dmatlist[[i]])<- ren.dupli(rownames(dmatlist[[i]]))
            
            return(list(d=dmatlist,tempProbs=tempprobs,conlimit=conlim,prob=prob,nhap=nhaps,rowindex=grouplist))
        }
    } else {
        
        
        tempprobs<-numeric(0)
        conlim<-Inf
        dmat<-clustersG(d,nhap)
        dmat<-geodist.corr(dmat)
        nm<-paste("net",1,sep="")
        grouplist<-list(1:nhap)
        dmatlist<-list(dmat)
        names(grouplist)<-nm
        names(dmatlist)<-nm
        
        numnet<-length(dmatlist)
        for(i in 1:numnet) rownames(dmatlist[[i]])<- ren.dupli(rownames(dmatlist[[i]]))
        
        
        return(list(d=dmatlist,tempProbs=tempprobs,conlimit=conlim,prob=0,nhap=nhap,rowindex=grouplist))
        
    }
    
}






#class Parsimnet

setClass(Class="Parsimnet", representation=representation(d="list",tempProbs="numeric",conlimit="numeric",prob="numeric",nhap="numeric",rowindex="list"))



#Set initialize Class Haplotype
setMethod("initialize", "Parsimnet", function(.Object,d=list(),tempProbs=numeric(0),conlimit=numeric(0),prob=numeric(0), nhap=numeric(0),rowindex=list())
{
    .Object@d<-d
    .Object@tempProbs<-tempProbs
    .Object@conlimit<-conlimit
    .Object@prob<-prob
    .Object@nhap<-nhap
    .Object@rowindex<-rowindex
    .Object
})



#fixing a bug
#Show method for Dna objects
setMethod("show","Parsimnet", function(object)
{
    
    cat("*** S4 Object of Class Parsimnet ***\n\n")
    cat("\nNumber of networks: \n")
    cat(length(object@d),"\n")
    cat("\nCalculated maximum connection steps at:",100*object@prob,"%:","\n")
    cat(object@conlimit,"\n")
    cat("\nNumber of haplotypes in each network: \n")
    cat(object@nhap,"\n")
    cat("\nNumber of intermediates in each network: \n")
    if(length(object@d)) cat(sapply(object@d,nrow)-object@nhap,"\n")
    cat("\nScore of each network (number of vertices minus one): \n")
    if(length(object@d)) cat(length(object),"\n")
    cat("\n\nslots of an object Parsimnet:\n", slotNames(object))
    
    
})


#new method 2016, cnrakt
#ren.dupli
ren.dupli<-function(x)
{
    
    rn<-x
    tb<-table(rn)
    if(any(tb>1))
    {
        dupl<-which(tb>1)
        
        for(t in 1:length(dupl))
        {
            namdup<-names(tb[dupl[t]])
            
            rn[rn==namdup][-1]<-paste(namdup,"_",2:tb[dupl[t]],sep="")
            
        }
        
        #warning("Duplicate names are renamed.")
    }
    
    rn
}


#re-run the codes after fixing bugs/update
#Generic parsimnet

setGeneric (
name= "parsimnet",
def=function(x,...)standardGeneric("parsimnet")
)


#parsimnet method for Dna objects
#Notes: ... additional arguments to phangorn's

setMethod(f="parsimnet", signature= "Dna", definition=function(x,indels="sic",prob=.95)
{
    
    

    if(length(prob)>1||(prob<0.01|prob>0.99)&&(!is.null(prob))) stop("probability must be NULL or a numeric vector of length 1, in the range [0.01,0.99]")
    h<-haplotypes::haplotype(x,indels=indels)
    d<-h@d
    
    
    
        netw<-parsnet(d,x@seqlengths[1],prob=prob)
        
        return(new("Parsimnet",d=netw$d,tempProbs=netw$tempProbs,conlimit=netw$conlimit,prob=netw$prob, nhap=netw$nhap,rowindex=netw$rowindex))
        
    
    
}
)


#parsimnet method for matrix objects

setMethod(f="parsimnet", signature= "matrix", definition=function(x,seqlength,prob=.95)
{
    
    if(length(prob)>1||(prob<0.01|prob>0.99)&&(!is.null(prob))) stop("probability must be NULL or a numeric vector of length 1, in the range [0.01,0.99]")
    
    h<-haplotypes::haplotype(x)
    d<-h@d
    netw<-parsnet(d,seqlength=seqlength,prob=prob)
    new("Parsimnet",d=netw$d,tempProbs=netw$tempProbs,conlimit=netw$conlimit,prob=netw$prob,nhap=netw$nhap,rowindex=netw$rowindex)
    
}
)



#parsimnet method for dist objects

setMethod(f="parsimnet", signature= "dist", definition=function(x,seqlength,prob=.95)
{
    if(length(prob)>1||(prob<0.01|prob>0.99)&&(!is.null(prob))) stop("probability must be NULL or a numeric vector of length 1, in the range [0.01,0.99]")
    h<-haplotypes::haplotype(x)
    d<-h@d
    netw<-parsnet(d,seqlength=seqlength,prob=prob)
    new("Parsimnet",d=netw$d,tempProbs=netw$tempProbs,conlimit=netw$conlimit,prob=netw$prob,nhap=netw$nhap,rowindex=netw$rowindex)
}
)




#plot method for Parsimnet objects

setMethod(f="plot", signature(x = "Parsimnet", y = "missing"), definition=function(x,y,net=1,inter.labels=FALSE,...)
{
    
    dots <- list(...)
    arg<-names(dots)
    netnum<-length(x@d)
    i<-net
    netd<-x@d
    nm<-names(netd)[i]
    net<-netd[[i]]
    nhap<-x@nhap[i]
    rn<-rownames(net)
    net[net!=1]<-0
    if(nrow(net)==1) net[is.na(net)]<-1
    g<-network(net)
    if(!inter.labels&&(length(rn)>nhap)) rn[(nhap+1):length(rn)]<-NA
    if(!any(arg=="label")) dots$label<-rn
    if(!any(arg=="usearrows")) dots$usearrows <- FALSE
    if(nrow(net)>1) { if(!any(arg=="mode")) dots$mode<-"kamadakawai" } else   if(!any(arg=="mode")) dots$mode<-"circle"
    if(!any(arg=="pad")) dots$pad <- 1
    if(!any(arg=="label.cex")) dots$label.cex<-0.75
    if(!any(arg=="vertex.cex")) dots$vertex.cex <- c(rep(0.8,nhap),rep(0.5,nrow(net)-nhap))
    if(!any(arg=="vertex.col"))dots$vertex.col<-c(rep(2,nhap),rep(4,nrow(net)-nhap))
    if(!any(arg=="main")) dots$main<-nm
    args <- c(list(x = g),dots)
    do.call(plot, args)
    
}
)

#bug in col, fixed, may 2021
#new function: pieplot
#Generic pieplot

setGeneric (
name= "pieplot",
def=function(x,y,...)standardGeneric("pieplot")
)



setMethod(f="pieplot", signature= c("Parsimnet","Haplotype"), definition=function(x,y,net=1,factors, coord = NULL,inter.labels=FALSE,interactive=FALSE, rex=1,...)
{
    
    #Some internal structures of pieplot is taken from package network, author Carter T. Butts.
    
    p<-x
    h<-y
    dots <- list(...)
    arg<-names(dots)
    
    ns<-sum(h@freq)
    g<-grouping(h,factors)$hapmat
    lev<-ncol(g)
    hi<-p@rowindex[[net]]
    g<-g[hi,]
    rmc<-(apply(g,2,sum)>0)
    g<-g[, rmc]
    
    nhap<-p@nhap[net]
    
    d<- p@d[[net]]
    rn<-rownames(d)
    rownames(d)<-rn
    colnames(d)<-rn
    nm<-names(p@d[net])
    
    d[d!=1]<-0
    if(nrow(d)==1) net[is.na(d)]<-1
    
    d[lower.tri(d)]<-0
    con<-which(d==1,arr.ind=TRUE)
    
    nw<-network::network(d)
    
    if(!any(arg=="col")) dots$col<- rainbow(lev)
    dots$col<- dots$col[rmc]
    if(!inter.labels&&(length(rn)>nhap)) rn[(nhap+1):length(rn)]<-NA
    if(!any(arg=="label")) dots$label<-rn
    if(nrow(d)>1) { if(!any(arg=="mode")) dots$mode<-"fruchtermanreingold" } else   if(!any(arg=="mode")) dots$mode<-"circle"
    if(!any(arg=="pad")) dots$pad <- 0.5
    if(!any(arg=="label.cex")) dots$label.cex<-1
    if(!any(arg=="edge.col")) dots$edge.col<-1
    if(!any(arg=="edge.lwd")) dots$edge.lwd<-1
    if(!any(arg=="edge.lty")) dots$edge.lty<-1
    if(!any(arg=="vertex.col")) dots$vertex.col<-1
    if(!any(arg=="label.col")) dots$label.col<-1
    if(!any(arg=="label.pos")) dots$label.pos<-0
    if(!any(arg=="label.pad")) dots$label.pad<-1
    if(!any(arg=="displaylabels")) dots$displaylabels <-TRUE
    if(!any(arg=="xlab")) dots$xlab = ""
    if(!any(arg=="ylab")) dots$ylab = ""
    if(!any(arg=="vertex.sides")) dots$vertex.sides<-50
    if(!any(arg=="edges")) dots$edges<-200
    
    if(!any(arg=="radius"))
    {
        rad<-apply(g,1,sum)*rep(0.8,nhap)
        dots$radius<- rex*rad/max(rad)
    }
    if(!any(arg=="vertex.cex"))
    {
        
        if(nrow(d)-nhap) dots$vertex.cex <- c(rep(0.5,nrow(d)-nhap))*min(dots$radius)
    }
    
    dots$radius<-c(dots$radius,dots$vertex.cex)
    
    
    if(is.null(coord))
    {
        if(dots$mode=="circle") coords<-network::network.layout.circle(nw=nw, layout.par=dots$layout.par)
        if(dots$mode=="fruchtermanreingold") coords<-network::network.layout.fruchtermanreingold(nw=nw, layout.par=dots$layout.par)
        if(dots$mode=="kamadakawai") coords<-network::network.layout.kamadakawai(nw=nw, layout.par=dots$layout.par)
    } else coords<-coord
    
    cx <- coords[, 1]
    cy <- coords[, 2]
    
    xlim<-range(cx)
    ylim<-range(cy)
    xlim<-xlim+c(-max(dots$radius)/2-dots$pad ,max(dots$radius)/2+dots$pad)
    ylim<-ylim+c(-max(dots$radius)/2-dots$pad ,max(dots$radius)/2+dots$pad)
    plot(0, 0, xlim = xlim, ylim = ylim, type = "n", xlab =  dots$xlab, ylab =  dots$ylab, asp = 1, axes = FALSE,main=nm)
    
    x1<-cx[con[,1]];x2<-cx[con[,2]]
    y1<-cy[con[,1]];y2<-cy[con[,2]]
    
    for(i in 1:length(x1)) lines(c(x1[i],x2[i]), c(y1[i],y2[i]),col=dots$edge.col,lwd=dots$edge.lwd, lty= dots$edge.lty)
    
    
    par.args<-names(par())
    plot.network.default.args<-names(formals(network::plot.network.default))
    args1<-union(par.args, plot.network.default.args)
    
    floating.pie.args<-names(formals(plotrix::floating.pie))
    polygon.args<-names(formals(graphics::polygon))
    args2<-union(floating.pie.args, polygon.args)
    
    arg1_2<-setdiff(args1, args2)
    #localfloating.pie<-function()
    
    localfloating.pie<-function(bos) bos
    
    eval(parse(text=paste("localfloating.pie<-","function", "(","...,", paste( arg1_2,collapse=",",sep=" "),")","{","plotrix::floating.pie(...)","}",collapse="")))
    
    rds<-dots$radius
    cls<-dots$col
    dots2<-dots[is.na(match(names(dots),c("radius","col")))]
    
    for(i in 1: nhap)
    {
        wg<-g[i,]>0
        args.flp <- c(list(xpos = cx[i],ypos=cy[i],x=g[i,][wg], radius=rds[i],col= cls[wg],startpos=0),c(dots2,list(main=nm)))
        do.call(localfloating.pie,  args.flp)
    }
    
    
    if(nrow(d)-nhap) network::network.vertex(cx[-c(1:nhap)], cy[-c(1:nhap)], radius=dots$vertex.cex , sides = dots$vertex.sides, border = 1, col = dots$vertex.col, lty = 1, rot = 0, lwd = 1)
    
    
    
    
    if (dots$displaylabels & (!all(dots$label == ""))) {
        
        if (dots$label.pos == 0) {
            xhat <- yhat <- rhat <- rep(0, nrow(d))
            xoff <- cx - mean(cx)
            yoff <- cy - mean(cy)
            roff <- sqrt(xoff^2 + yoff^2)
            for (i in (1:nrow(d))) {
                ij <- unique(c(d[d[, 2] == i & d[, 1] != i, 1],
                d[d[, 1] == i & d[, 2] != i, 2]))
                ij.n <- length(ij)
                if (ij.n > 0) {
                    for (j in ij) {
                        dx <- cx[i] - cx[j]
                        dy <- cy[i] - cy[j]
                        dr <- sqrt(dx^2 + dy^2)
                        xhat[i] <- xhat[i] + dx/dr
                        yhat[i] <- yhat[i] + dy/dr
                    }
                    xhat[i] <- xhat[i]/ij.n
                    yhat[i] <- yhat[i]/ij.n
                    rhat[i] <- sqrt(xhat[i]^2 + yhat[i]^2)
                    if (!is.nan(rhat[i]) && rhat[i] != 0) {
                        xhat[i] <- xhat[i]/rhat[i]
                        yhat[i] <- yhat[i]/rhat[i]
                    }
                    else {
                        xhat[i] <- xoff[i]/roff[i]
                        yhat[i] <- yoff[i]/roff[i]
                    }
                }
                else {
                    xhat[i] <- xoff[i]/roff[i]
                    yhat[i] <- yoff[i]/roff[i]
                }
                if (is.nan(xhat[i]) || xhat[i] == 0)
                xhat[i] <- 0.01
                if (is.nan(yhat[i]) || yhat[i] == 0)
                yhat[i] <- 0.01
            }
            
            xhat <- xhat
            yhat <- yhat
        }
        else if (dots$label.pos < 5) {
            xhat <- switch(dots$label.pos, 0, -1, 0, 1)
            yhat <- switch(dots$label.pos, -1, 0, 1, 0)
        }
        else if (dots$label.pos == 6) {
            xoff <- cx - mean(cx)
            yoff <- cy- mean(cy)
            roff <- sqrt(xoff^2 + yoff^2)
            xhat <- xoff/roff
            yhat <- yoff/roff
        }
        else {
            xhat <- 0
            yhat <- 0
        }
        os <- par()$cxy * mean(dots$label.cex, na.rm = TRUE)
        lw <- strwidth(dots$label, cex = dots$label.cex)/2
        lh <- strheight(dots$label, cex = dots$label.cex)/2
        
        text(cx + xhat * rds + (lh * dots$label.pad + lw) * ((xhat > 0) - (xhat < 0)), cy + yhat * rds + (lh * dots$label.pad + lh) * ((yhat > 0) - (yhat < 0)), dots$label, cex = dots$label.cex, col = dots$label.col, offset = 0)
    }
    
    
    
    if(interactive)
    
    {
        "%iin%" <- function(x, int) (x >= int[1]) & (x <= int[2])
        if (interactive && ((length(cx) > 0))) {
            os <- c(0.2, 0.4) * par()$cxy
            textloc <- c(min(cx) - dots$pad, max(cy) + dots$pad )
            tm <- "Select a vertex to move, or click \"Finished\" to end."
            tmh <- strheight(tm)
            tmw <- strwidth(tm)
            text(textloc[1], textloc[2], tm, adj = c(0, 0.5))
            fm <- "Finished"
            finx <- c(textloc[1], textloc[1] + strwidth(fm))
            finy <- c(textloc[2] - 3 * tmh - strheight(fm)/2, textloc[2] -
            3 * tmh + strheight(fm)/2)
            finbx <- finx + c(-os[1], os[1])
            finby <- finy + c(-os[2], os[2])
            rect(finbx[1], finby[1], finbx[2], finby[2], col = "white")
            text(finx[1], mean(finy), fm, adj = c(0, 0.5))
            clickpos <- unlist(locator(1))
            if ((clickpos[1] %iin% finbx) && (clickpos[2] %iin% finby)) {
                cl <- match.call()
                cl$interactive <- FALSE
                cl$coord <- cbind(cx, cy)
                return(eval.parent(cl))
            }
            else {
                clickdis <- sqrt((clickpos[1] - cx)^2 + (clickpos[2] -
                cy)^2)
                selvert <- match(min(clickdis), clickdis)
                if (all(dots$label == ""))
                dots$label <- 1:nrow(d)
                rect(textloc[1], textloc[2] - tmh/2, textloc[1] +
                tmw, textloc[2] + tmh/2, border = "white", col = "white")
                tm <- "Where should I move this vertex?"
                tmh <- strheight(tm)
                tmw <- strwidth(tm)
                text(textloc[1], textloc[2], tm, adj = c(0, 0.5))
                fm <- paste("Vertex", dots$label[selvert], "selected")
                finx <- c(textloc[1], textloc[1] + strwidth(fm))
                finy <- c(textloc[2] - 3 * tmh - strheight(fm)/2,
                textloc[2] - 3 * tmh + strheight(fm)/2)
                finbx <- finx + c(-os[1], os[1])
                finby <- finy + c(-os[2], os[2])
                rect(finbx[1], finby[1], finbx[2], finby[2], col = "white")
                text(finx[1], mean(finy), fm, adj = c(0, 0.5))
                clickpos <- unlist(locator(1))
                cx[selvert] <- clickpos[1]
                cy[selvert] <- clickpos[2]
                cl <- match.call()
                cl$coord <- cbind(cx, cy)
                return(eval.parent(cl))
                
                
            }
        }
        
        
    }
    invisible(cbind(cx, cy))
}
)



#new function: pielegend
#Generic pielegend

setGeneric (
name= "pielegend",
def=function(p,h,...)standardGeneric("pielegend")
)


setMethod(f="pielegend", signature= c("Parsimnet","Haplotype"), definition=function(p,h,net=1,factors,...)
{
    dots <- list(...)
    arg<-names(dots)
    ns<-sum(h@freq)
    g<-grouping(h,factors)$hapmat
    lev<-ncol(g)
    hi<-p@rowindex[[net]]
    g<-g[hi,]
    rmc<-(apply(g,2,sum)>0)
    fac<-levels(factor(factors))[rmc]
    
    if(!any(arg=="col")) dots$col<- rainbow(lev)
    dots$col<- dots$col[rmc]
    
    if(!any(arg=="fill")) dots$fill <-dots$col
    if(!any(arg=="legend")) dots$legend<-fac
    if(!any(arg=="x")&&!any(arg=="y")) dots$x<-"topright"
    do.call(graphics::legend, dots)
    
})



#Update
#length method for Parsimnet objects

setMethod(f="length", signature= "Parsimnet", definition=function(x)
{
    
    sapply(x@d,nrow)-1
    
}
)

#new method September 2016, cnrakt
#names method for Parsimnet objects

setMethod("names","Parsimnet", function(x)
{
    names(x@d)
})

#new method September 2016, cnrakt
#names replace method for Dna objects

setReplaceMethod("names","Parsimnet", function(x,value)
{
    if(is.numeric(value)) value<-as.character(value)
    names(x@d)<-value
    x
})
#new method September 2016, cnrakt
#rownames method for Parsimnet objects

setMethod("rownames","Parsimnet", function(x)
{
    n<-length(x@d)
    rownameslist<-vector("list",n)
    for(i in 1:n)
    {
        rownameslist[[i]]<-rownames(x@d[[i]])
        
    }
    rownameslist
})

#new method September 2016, cnrakt
#rownames replace method for Dna objects

setReplaceMethod("rownames","Parsimnet", function(x,value)
{
    
    n<-length(x@d)
    for(i in 1:n)
    {
        rn<-as.character(value[[i]])
        rownames(x@d[[i]])<-rn
        colnames(x@d[[i]])<-rn
    }
    
    return(x)
})


#Coerce Parsimnet object to a list

setMethod("as.list","Parsimnet", function(x)
{
    l<-list(x@d,x@tempProbs, x@conlimit,x@prob,x@nhap,x@rowindex)
    names(l)<-c("d","tempProbs","conlimit","prob","nhap","rowindex")
    l
})



#as.network method for parsimnet objects
#setOldClass("as.network")
setMethod(f="as.network", signature= c("Parsimnet"), definition=function(x,net=1,...)
{
    
    net<-x@d[[net]]
    rn<-rownames(net)
    nhap<-x@nhap[net]
    net[net!=1]<-0
    if(nrow(net)==1) net[is.na(net)]<-1
    g<-network(net,...)
    g
})




#as.networx
#setOldClass("as.networx")
setMethod(f="as.networx", signature= c("Parsimnet"), definition=function(x,net=1,...)
{
    
    nhap<-x@nhap[net]
    d<-x@d[[net]][1:nhap,1:nhap]
    rn<-rownames(d)
    n<-NJ(d)
    spp<-as.splits(n,...)
    netp<-as.networx(spp)
    netp
})



#new function
#internal function geodist.corr (for fixing a bug in parsimnet)
geodist.corr<-function(dmat,...)
{
    net<-dmat
    net[net!=1]<-0
    if(nrow(net)==1) net[is.na(net)]<-1
    g<-geodist(network(net),...)$gdist
    rownames(g)<-rownames(dmat)
    colnames(g)<-colnames(dmat)
    g
}



#new internal function
#internal function has.cycle
has.cycle<-function(d)
{
    
    d<-as.matrix(d)
    dupd<-!duplicated(d)
    d<-d[dupd,dupd]
    d[d!=1]<-0
    d[lower.tri(d)]<-0
    return(1+sum(d)-nrow(d))
    
}



#new internal function
#internal function:  ext.nodes
ext.nodes<-function(d)
{
    d<-as.matrix(d)
    dupd<-!duplicated(d)
    if(sum(dupd)<3) return(rownames(d))
    d<-d[dupd,dupd]
    d[d!=1]<-0
    lab<-(rownames(d))
    r.s<-rowSums(d)
    extseqslabs<-match(lab[r.s==1],rownames(d))
    return(rownames(d)[extseqslabs])
}




#new function: homopoly
#Generic homopoly

setGeneric (
name= "homopoly",
def=function(x,...)standardGeneric("homopoly")
)



setMethod(f="homopoly", signature= c("Dna"), definition=function(x,indels="sic",...)
{
    
    indmet<- c("sic","5th","missing")
    matchmet <- pmatch(indels, indmet)
    if (is.na(matchmet))
    stop("invalid indel coding method", paste("", indels))
    indels<-indmet[matchmet]
    x<-haplotypes::as.dna(haplotypes::haplotype(x,indels=indels))
    dat<-haplotypes::as.phyDat(x,indels=indels)
    p<-parsimnet(x,indels=indels,prob=NULL,...)
    tree<-ape::nj(p@d[[1]][1:p@nhap,1:p@nhap])
    ci<-phangorn::CI(tree,dat,sitewise=TRUE)
    homop<-which(ci<1)
    
    ind<-homop[homop>ncol(x)]
    
    if(length(ind))
    {
        indn<-ind-ncol(x)
        ind<-ci[ind]
        names(ind)<-indn
    }
    
    su<-homop[homop<=ncol(x)]
    sun<-su
    su<-ci[su]
    names(su)<-sun
    
    
    
    return(c(indels=list(ind), subs=list(su)))
    
    
})







