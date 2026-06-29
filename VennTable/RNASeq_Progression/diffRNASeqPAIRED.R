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

##
## Then get paired samples

progPaired=prog[prog$FIRST.RNASEQ=="X" & prog$NORMAL.RNASEQ=="X",]
pairedSamples=c(progPaired$NORMAL,progPaired$FIRST)
nPairs=nrow(progPaired)

ss.paired=intersect(pairedSamples,colnames(ds))
ds=ds[,ss.paired]

y <- DGEList(counts=ds)

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

#plotMDS(y)

pNF=data.frame(Patient=seq(nPairs),Sample=progPaired[,1],Tissue=rep("NF",4))
pWD=data.frame(Patient=seq(nPairs),Sample=progPaired[,2],Tissue=rep("WD",4))
key=rbind(pNF,pWD)

Patient=factor(key$Patient)
Tissue=factor(key$Tissue)

design <- model.matrix(~Patient+Tissue)

y <- estimateDisp(y, design, robust=T)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)

aa=topTags(lrt,n=nrow(ds),adjust="fdr")

qCut=1.01
i.sig=which(aa$table$FDR<qCut)
sig.genes=rownames(aa$table)[i.sig]

rescale=mean(log2(colSums(ds[,pairedSamples])))-mean(log2(colSums(cpm(y))))
ldn=(log2(cpm(y))+rescale)
dn=2^ldn

meanSig=2^t(apply(log2(dn[sig.genes,]),1,function(x){tapply(x,Tissue,mean)}))

lfc=aa$table$logFC
FC=ifelse(lfc<0,-(2^-lfc),2^lfc)
CPM=2^aa$table$logCPM

ans=data.frame(
    FDR=aa$table$FDR[i.sig],
    FC=FC[i.sig],
    dn[sig.genes,]
    )

write.xls(ans,file=cc("rnaSEQ_FIRST_WD_vs_NF_Paired__FDR__v2",qCut,".txt"))

