data(cghUBM)
data(progressionSet)
data(sampleTable)

threshold=0.3

progM=progressionSet[progressionSet$PROG!="",]

integrateSamples=sampleTable$BIO_ID[sampleTable$aCGH_U133A=="Y"]

#
# 2024-01-30
# Get sample manifest
#

manifest=sampleTable[sampleTable$aCGH_U133A=="Y",]
typeCount=tibble::tibble(manifest) %>% dplyr::count(TYPE)
openxlsx::write.xlsx(list(Types=typeCount,Manifest=manifest),"progressionBlocks_WDvsDD_Expr_v1__Manifest.xlsx")


amplifed=(ubm$a0>threshold & ubm$d0<threshold)
deleted=(ubm$a0<threshold & ubm$d0>threshold)

amplifedPairs=(!amplifed[,progM$FIRST] & amplifed[,progM$PROG])
deletedPairs=(!deleted[,progM$FIRST] & deleted[,progM$PROG])

numAmpPairs=apply(amplifedPairs,1,sum)
numDelPairs=apply(deletedPairs,1,sum)

jAmp=which(numAmpPairs>0)
jDel=which(numDelPairs>0)

collapsePairs<-function(x){
  cc(progM$FIRST[x],progM$PROG[x],collapse=";")
}

pairsAmp=apply(amplifedPairs,1,collapsePairs)
pairsDel=apply(deletedPairs,1,collapsePairs)

Mamp=ifelse(amplifed[jAmp,integrateSamples],"A","")
Mdel=ifelse(deleted[jDel,integrateSamples],"D","")

type=factor(sampleTable$TYPE[sampleTable$aCGH_U133A=="Y"],levels=c("WD","DD"))


M1=data.frame(
    BreakNo=jAmp,
    ubm$breaks[jAmp,1:3],
    Event="Amp",
    NumPaired=numAmpPairs[jAmp],
    Pairs=pairsAmp[jAmp])

M2=data.frame(
    BreakNo=jDel,
    ubm$breaks[jDel,1:3],
    Event="Del",
    NumPaired=numDelPairs[jDel],
    Pairs=pairsDel[jDel])
M12=rbind(M1,M2)

Mtbl=rbind(Mamp,Mdel)

M12$CGH.FT.pVal=apply(Mtbl,1,function(x){fisher.test(table(x,type))$p.value})
M12$CGH.FT.FDR="NA"
M12$CGH.frac=t(apply(Mtbl,1,function(x){tapply(x!="",type,mean,na.rm=T)}))


MM=cbind(M12,Mtbl)
ii=order(M12$Chr,M12$Start)
MM=MM[ii,]

Mtbl=Mtbl[ii,]

eventSig=apply(Mtbl,1,paste0,collapse=",")

MM.o=MM

dups=which(duplicated(eventSig))
for(ii in dups) {

    print(ii)

    MM$BreakNo[ii]=paste(MM$BreakNo[ii-1],MM$BreakNo[ii],sep=",")
    MM$Start[ii]=MM$Start[ii-1]
    MM$BreakNo[ii-1]=""
    MM$Start[ii-1]=-1

}

MMdd=MM[MM$BreakNo!="",]
MMdd$CGH.FT.FDR=p.adjust(MMdd$CGH.FT.pVal,"fdr")

qCut=0.05
MMsig=MMdd[MMdd$CGH.FT.FDR<qCut,]

data(u133a)
ds=u133a$ds[,colnames(Mtbl)]
require(GenomicRanges)
u133a$probesRanges=GRanges(
    seqnames=u133a$annoteHG18$CHR,
    IRanges(
        start=u133a$annoteHG18$START,
        end=u133a$annoteHG18$END
        ))

tbl=NULL
for(ii in seq(nrow(MMsig))) {

    chrom=paste0("chr",MMsig$Chr[ii])
    start=MMsig$Start[ii]
    end=MMsig$End[ii]
    event=MMsig$Event[ii]

    blockRange=GRanges(seqnames=chrom,IRanges(start=start,end=end))

    ff=findOverlaps(blockRange,u133a$probesRanges)
    overlap=pintersect(blockRange[queryHits(ff)],u133a$probesRanges[subjectHits(ff)])
    pct.overlap=width(overlap)/width(u133a$probesRanges[subjectHits(ff)])
    pij=pct.overlap>0.75

    probes=rownames(u133a$annoteHG18)[subjectHits(ff)[pij]]

    if(len(probes)==0)
        next

    for(pp in seq(len(probes))) {
        cat(ii,chrom,start,end,event,pp,probes[pp],"\n")

        dat=data.frame(
            TYPE=type,
            Expr=ds[probes[pp],]
            )

        xx=summary(aov(Expr ~ TYPE,dat=dat))
        pValue=xx[[1]][1,5]
        mSig=tapply(dat$Expr,dat$TYPE,mean)

        tbii=list(
            BLOCK=ii,
            MMsig[ii,1:10],
            GENE=u133a$annote[probes[pp],1],
            Probe=probes[pp],
            Pvalue=pValue,
            log2FC=diff(mSig)
            )

        if(is.data.frame(tbl)){
            tbl=rbind(tbl,data.frame(tbii))
        } else {
            tbl=data.frame(tbii)
        }

    }

}

# Only events with consistent sign
jj=(tbl$log2FC>0 & tbl$Event=="Amp") | (tbl$log2FC<0 & tbl$Event=="Del")
tbl=tbl[jj,]

# Only events with real genes
tbl=tbl[!is.na(tbl$GENE),]

tbl$FDR=p.adjust(tbl$Pvalue,"fdr")
tbl$FC=ifelse(tbl$log2FC<0,-2^(-tbl$log2FC),2^tbl$log2FC)
tbl=tbl[order(tbl$Pvalue),]

VERSION="v1"
OTAG=cc("progressionBlocks_WDvsDD_Expr",VERSION)

library(openxlsx)
write.xlsx(
    list(chr13=tbl[tbl$Chr==13,], chr8=tbl[tbl$Chr==8,], ALL=tbl),
    cc(OTAG,".xlsx"),
    rowNames=FALSE
)

pdf(file=cc("boxplots",OTAG,".pdf"),height=11,width=8.5)
par(mfrow=c(3,2))


stab=list()

for(jj in seq(nrow(tbl))) {
    cat(jj,"\n")
    probes_pp=tbl$Probe[jj]
    dat=data.frame(
            Type=type,
            Expr=ds[probes_pp,]
            )

    minStats=min(sapply(tapply(dat$Expr,dat$Type,boxplot.stats),function(x){min(x$stats)}))
    minOut=min(sapply(tapply(dat$Expr,dat$Type,boxplot.stats),function(x){min(x$out)}))
    yMin=min(minOut,minStats)

    maxStats=max(sapply(tapply(dat$Expr,dat$Type,boxplot.stats),function(x){max(x$stats)}))
    maxOut=max(sapply(tapply(dat$Expr,dat$Type,boxplot.stats),function(x){max(x$out)}))
    yMax=max(maxOut,maxStats)

    BLOCK=paste0("B",tbl$BLOCK[jj])
    if(is.null(stab[[BLOCK]])){
        stab[[BLOCK]]=dat[,1,drop=F]
    }
    stab[[BLOCK]]=cbind(stab[[BLOCK]],dat[,2,drop=F])
    colnames(stab[[BLOCK]])[ncol(stab[[BLOCK]])]=paste0(u133a$annote[probes_pp,1],"|",probes_pp)

    bb=boxplot(Expr ~ Type, dat=dat, outline=F,
        ylim=c(yMin,yMax),
        main=paste(u133a$annote[probes_pp,1],"|",probes_pp,"\n",
            paste0("Chr",tbl$Chr[jj],":",tbl$Start[jj],"-",tbl$End[jj])),
        sub=paste("FDR =",prettyNum(tbl$FDR[jj]),"FC =",prettyNum(tbl$FC[jj])))

    xj=jitter(as.numeric(dat$Type))
    points(xj,dat$Expr,pch=19,col=8)

}


dev.off()




