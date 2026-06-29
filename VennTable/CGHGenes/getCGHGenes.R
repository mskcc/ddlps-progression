data(cghGeneMatrix)
cgh_probes=rownames(geneCGH)

dd=read.delim("../U133A_NFvsWD/u133A_WDrs_vs_NF_Pfilter_ALL_v4.txt")
dd$RAE.GENE=NA

for(ii in seq(nrow(dd))){

    gene=dd$SYMBOL[ii]
    jj=grep(paste0(":",gene,"$"),cgh_probes)
    cat("Gene,#matches =",gene,len(jj),"\n")

    if(len(jj)>0){
        evtSign=sign(dd$FC[ii])
        cghCalls=sign(as.numeric(as.matrix(geneCGH[jj,])))
        cghCalls=cghCalls[!is.na(cghCalls)]
        correctSign=sum(cghCalls==evtSign)/len(cghCalls)
        discordantSign=sum(cghCalls==(-1)*(evtSign))/len(cghCalls)
        if(correctSign > discordantSign){
            dd$RAE.GENE[ii]=round(correctSign,3)
        } else {
            dd$RAE.GENE[ii]=round(-1,0)
        }
    }
}

# Get rid of multi genes
dd=dd[-grep(";",dd$SYMBOL),]

# Find 12q genes

require(IRanges)
data(u133a)

ii.12=which(u133a$annoteHG18[dd$X,1]=="chr12")
u133a.12=IRanges(start=u133a$annoteHG18[dd$X,][ii.12,]$START,end=u133a$annoteHG18[dd$X,][ii.12,]$END)
data(chr12qEvent)
chr12q=IRanges(start=chr12qEvent$Start,end=chr12qEvent$End)
ii.12q=which(u133a.12 %over% chr12q)


dd$Chr12q=""
dd$Chr12q[ii.12][ii.12q]="X"

ans=data.frame(SYMBOL=dd$SYMBOL,RAE.GENE=dd$RAE.GENE,Chr12q=dd$Chr12q)
write.xls(ans,file="geneCGHvsU133aConcordence_v3.txt")
