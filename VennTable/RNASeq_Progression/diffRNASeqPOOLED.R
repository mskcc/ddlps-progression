require(edgeR)
require(limma)

data(rnaSeq)
data(sampleTable)
ds=rnaSeq$ds

##
## Get just progression FIRST Samples
##

data(progressionSet)
prog=progressionSet

ss.first=intersect(c(prog$FIRST,prog$NORMAL),colnames(ds))

ds=ds[,ss.first]
subTypes=sampleTable$TYPE[match(colnames(ds),sampleTable$BIO_ID)]

group=factor(subTypes,levels=c("NF","WD"))

y <- DGEList(counts=ds,group=group)

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]


y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
et <- exactTest(y)

aa=topTags(et,n=nrow(ds),adjust="fdr")

#######
#
# Continue
#



qCut=1.01
i.sig=which(aa$table$FDR<1.01)
sig.genes=rownames(aa$table)[i.sig]

rescale=mean(log2(colSums(ds[,ss.first])))-mean(log2(colSums(cpm(y))))
dn=2^(log2(cpm(y))+rescale)

meanSig=2^t(apply(log2(dn[sig.genes,]),1,function(x){tapply(x,group,mean)}))

lfc=aa$table$logFC
FC=ifelse(lfc<0,-(2^-lfc),2^lfc)

ans=data.frame(
    FDR=aa$table$FDR[i.sig],
    FC=FC[i.sig],
    meanSig
    )


write.xls(ans,file=cc("rnaSEQ_FIRST_WD_vs_NF_Pop__FDR",qCut,".txt"))

