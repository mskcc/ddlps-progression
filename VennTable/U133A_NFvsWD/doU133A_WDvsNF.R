library(limma)
library(gdata)

data(sampleTable)
data(u133a)

bioids=sampleTable$BIO_ID[sampleTable$U133A=="Y" & sampleTable$TYPE %in% c("WD","NF")]
ds=u133a$ds[,bioids]
dp=u133a$dp[,bioids]
type=factor(sampleTable[colnames(ds),"TYPE"])

filter.P <- function(dp,cl) {
  cutP = table(cl)/4
  ig=apply(dp,1,P.filter,cl,cutP)
  return(ig)
}

P.filter <- function(calls,fact,cutP) {
  numP = tapply(calls=="P",fact,sum)
  return(any(numP>cutP))
}



gg=filter.P(dp,type)
cat(c(">>>Num.P.genes=",sum(gg),"\n"))

ds=ds[gg,]
dp=dp[gg,]

design=model.matrix(~0+type)
colnames(design)=gsub("type","",colnames(design))
fit=lmFit(ds,design)

col.design=colnames(design)

N.t=length(col.design)
cm=matrix(0,nrow=N.t,ncol=N.t*(N.t-1)/2)

dimnames(cm)=list(
    Levels=colnames(design),
    Contrasts=seq(ncol(cm)))

##
## All pairwise contrasts
##

con=1
for(i in 1:(N.t-1)) {
  for(j in (i+1):N.t) {
    cm[col.design[i],con]=-1
    cm[col.design[j],con]=1
    colnames(cm)[con]=paste(col.design[j],col.design[i],sep="-")
    con = con + 1
  }
}


fit2=contrasts.fit(fit,cm)
fit2=eBayes(fit2)

fit2$ds=ds

aa=topTable(fit2,n=nrow(ds))

ig=nrow(ds)

ans=aa[1:ig,]
ans$FC=ifelse(ans$logFC<0,-2^(-ans$logFC),2^ans$logFC)

ans=data.frame(ans,u133a$annote[rownames(ans),])

ans=ans[,c("SYMBOL","adj.P.Val","FC","CHR","CHRLOC")]
gi=rownames(ans)
meanSig=data.frame(mean=2^t(apply(ds[gi,],1,function(x){tapply(x,type,mean)})))
ans=cbind(ans,meanSig)

write.xls(ans,file=cc("u133A_WDrs_vs_NF_Pfilter_ALL_v4.txt"))
