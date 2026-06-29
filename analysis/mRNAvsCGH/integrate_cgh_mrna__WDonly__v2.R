data(u133a)
data(cghEventsPerBlock)
data(sampleTable)

bioId=sampleTable$BIO_ID[sampleTable$aCGH_U133A=="Y"]

ds=u133a$ds[,bioId]
cghEvents=cghEvents[bioId,]

eventTbl=readr::read_csv("data/db/progressionEventTable.csv")
rois=c("1q21-24","6q23-25","19p13.3","19p13.2","19p13.13")
eventTbl=dplyr::filter(eventTbl,(Band %in% rois))

regions=eventTbl$Region

tbl=NA

for(ii in seq(len(regions))) {

    chrom=paste0("chr",strsplit(regions[ii],"\\.")[[1]][2])
    start=as.numeric(strsplit(regions[ii],"\\.")[[1]][3])
    end=as.numeric(strsplit(regions[ii],"\\.")[[1]][4])
    eventType=strsplit(regions[ii],"\\.")[[1]][1]

#
#   Only get probes in event area
#
    probes=unique(
        rownames(u133a$annoteHG18)[
        (chrom==u133a$annoteHG18$CHR & start<=u133a$annoteHG18$START & end>=u133a$annoteHG18$END)
        ])

#    probes=unique(rownames(u133a$annoteHG18))

    for(pp in seq(len(probes))) {
        cat(ii,regions[ii],pp,probes[pp],"\n")

        dat=data.frame(
            TYPE=sampleTable[bioId,"TYPE"],
            Expr=ds[probes[pp],],
            CopyNum=cghEvents[,regions[ii]])

        samplesToUse=(dat$TYPE=="WD")
        dat=dat[samplesToUse,]

        dat$Status=NA

        dat$Status[dat$TYPE=="WD" & dat$CopyNum=="X"]=paste0("WD_",eventType)
        dat$Status[dat$TYPE=="WD" & dat$CopyNum!="X"]="WD.0"

        dat$Status=factor(dat$Status,
            levels=c("WD.0",paste0("WD_",eventType)))

        xx=summary(aov(Expr ~ Status, dat=dat))
        pValue=xx[[1]][1,5]
        mSig=tapply(dat$Expr,dat$Status,mean)

        tbii=list(
            Region=regions[ii],
            Event=eventType,
            GENE=u133a$annote[probes[pp],1],
            Probe=probes[pp],
            Pvalue=pValue,
            log2FC=as.numeric((mSig[len(mSig)]-mSig[len(mSig)-1]))
            )

        if(is.data.frame(tbl)){
            tbl=rbind(tbl,data.frame(tbii))
        } else {
            tbl=data.frame(tbii)
        }
        #print(nrow(tbl))

    }
}


# Only events with consistent sign
#jj=(tbl$log2FC>0 & tbl$Event=="Amp") | (tbl$log2FC<0 & tbl$Event=="Del")
#tbl=tbl[jj,]

tbl$FDR=p.adjust(tbl$Pvalue,"fdr")
tbl=tbl[order(tbl$Pvalue),]
tbl$FC=ifelse(tbl$log2FC<0,-2^(-tbl$log2FC),2^tbl$log2FC)

#stop("88")

VERSION="v7"
OTAG=cc("mRNA_vs_CopyNumber_by_blocks___WDonly",VERSION)
save(tbl,file=cc("diff",OTAG,".rda"),compress=T)

require(tidyverse)

source("Rlib/annotation.R")
cq<-paste0
ann=getAnnoteTable("Rlib/hgu133a.sqlite")
ann=tibble(ann) %>% select(GENE=symbol,Description=gene_name) %>% distinct

tbl2=tibble(tbl) %>%
    left_join(eventTbl) %>%
    left_join(ann,relationship = "many-to-many") %>%
    distinct(Region,GENE,Probe,.keep_all=T) %>%
    select(Band,Event,Region,GENE,Description,everything()) %>%
    rename(`Gene symbol`=GENE) %>%
    arrange(`Gene symbol`)

cghGeneTbl=readxl::read_xlsx("../../tables/CGHEventTable/cghGeneTable___Manuscript___SelectRegions1.xlsx") %>% select(1:3)
bandEvent=read_csv("bandEvent.csv")

# Rewrite legacy HGNC aliases in cghGeneTbl to current symbols so the join
# matches tbl2's u133a-derived symbols. Without this, e.g. MAP3K7IP2 in
# cghGeneTbl never matches TAB2 in tbl2 and the row gets NA stats.
alias_map=AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys=unique(cghGeneTbl$geneSymbol),
        columns=c("SYMBOL"),
        keytype="ALIAS") %>%
    filter(!is.na(SYMBOL), ALIAS != SYMBOL) %>%
    distinct(ALIAS, .keep_all=TRUE) %>%
    rename(legacy=ALIAS, modern=SYMBOL)

cghGeneTbl=cghGeneTbl %>%
    left_join(alias_map, by=c(geneSymbol="legacy")) %>%
    mutate(geneSymbol=coalesce(modern, geneSymbol)) %>%
    select(-modern)

tbl3=cghGeneTbl %>%
    left_join(tbl2,by=c(Event="Region",geneSymbol="Gene symbol")) %>%
    select(-Band) %>%
    left_join(bandEvent) %>%
    select(`Gene Symbol`=geneSymbol,Description=description,Band,Probe,Pvalue,log2FC,FDR,FC) %>%
    group_split(Band)

SET_TAG="WDpEvent_vs_WD_without_v2"

names(tbl3)=map_vec(tbl3,~.$Band[1])
tbl3=tbl3 %>% map(~select(.,-Band)) %>% map(data.frame,check.names=F)
openxlsx::write.xlsx(tbl3,paste0("suppTable_4_with_U133A__",SET_TAG,".xlsx"))

#
# === SUPPLEMENTAL TABLE 7 (drop-in replacement) ===
# Supplemental Table 7 is the 6q23-25 band, restricted to genes overexpressed
# at fold change > 1.4 in WDLS with vs. without the amplicon, ordered by fold
# change. The published table's FC and P value are exact, but its FDR column was
# carried over from an earlier run with a different probe set and is not
# reproducible from the current pipeline. This block re-emits the SAME rows and
# columns with the now-reproducible BH FDR (computed by p.adjust over all probes
# in the five amplicon regions, line 82 -- n=778). It is written to a NEW file
# (suppTable7.xlsx) and does NOT touch the workbook written above.
#
suppTable7=tbl3[["6q23-25"]] %>%
    filter(!is.na(FC), FC > 1.4) %>%
    arrange(desc(FC)) %>%
    select(Gene=`Gene Symbol`,Probe,`Fold change`=FC,`P value`=Pvalue,FDR)

write_xlsx(suppTable7,"suppTable7.xlsx")
