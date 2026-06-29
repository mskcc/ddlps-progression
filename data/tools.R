require(stringr)
fixCGHSlideNames<-function(x){
    y=str_match(x,"(.*)(_CGH)")[,2]
    gsub("_S\\d+$","",y)
}

#cnv=read.delim("CNV.file",skip=1)
source("raeCNV.R")
require(IRanges)
getPolymorphicCNV<-function(breaks,cnv=raeCNV) {
    polymorphic=NULL
    for(ch in unique(breaks$Chr)) {

        idx=which(breaks$Chr==ch)
        brks=IRanges(start=breaks$Start[idx],end=breaks$End[idx])
        poly=IRanges(start=cnv$Start[cnv$Chr==ch],end=cnv$End[cnv$Chr==ch])
        cov1=lapply(1:length(brks),function(b) coverage(poly,shift=-start(brks[b,])+1,width=width(brks[b,])))
        cov2=unlist(lapply(cov1,function(cc) sum(width(slice(cc,lower=1,includeLower=TRUE)))/length(cc)))

        polymorphic=c(polymorphic,idx[cov2>0.5])

    }
    return(polymorphic)
}

imputePolymorphicSomatic<-function(polymorphic,breaks,xx) {

    som=(1:nrow(breaks))[-polymorphic]
    near=som[nearest(
        IRanges(start=polymorphic,end=polymorphic),
        IRanges(start=som,end=som)
    )]
    xx[polymorphic]=xx[near]
    return(xx)

}



