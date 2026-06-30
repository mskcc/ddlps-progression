source("R/helpers.R")

data(u133a)
data(cghEventsPerBlock)
data(sampleTable)

bioId=sampleTable$BIO_ID[sampleTable$aCGH_U133A=="Y"]

ds=u133a$ds[,bioId]
cghEvents=cghEvents[bioId,]

#
# Get gene symbol, name/desc annotation
#

geneAnn=AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        columns=c("SYMBOL","GENENAME"),
        keys=unname(unique(as.character(org.Hs.eg.db::org.Hs.egSYMBOL))),
        keytype="SYMBOL"
    ) |>
    dplyr::distinct(SYMBOL,.keep_all=T)

#
# Which comparision (subtype)
#

subType="DD"

SET_TAG=paste0(subType,"pEvent_vs_",subType,"_without")

TABLE_NO="12"
#
#
#


#
# Get ROI's for comparisoin
#

sTableFile="manuscriptTables/Supplemental Table 12.Genes encoded on DD associated CNA.xlsx"
sheets=readxl::excel_sheets(sTableFile)

eventTbl=purrr::map(sheets,~readxl::read_xlsx(sTableFile,sheet=.)) |>
    purrr::map(dplyr::distinct,Band,Event) |>
    dplyr::bind_rows() |>
    dplyr::rename(Region=Event)

regions=eventTbl$Region

#
#
#

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

        samplesToUse=(dat$TYPE==subType)
        dat=dat[samplesToUse,]

        dat$Status=NA

        dat$Status[dat$TYPE==subType & dat$CopyNum=="X"]=paste0(subType,"_",eventType)
        dat$Status[dat$TYPE==subType & dat$CopyNum!="X"]=paste0(subType,".0")

        dat$Status=factor(
                            dat$Status,
                            levels=c(
                                paste0(subType,".0"),
                                paste0(subType,"_",eventType)
                            )
                        )

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

VERSION="v8"
OTAG=paste0("mRNA_vs_CopyNumber_by_blocks___",subType,"only_",VERSION)
save(tbl,file=cc("diff",OTAG,".rda"),compress=T)

require(tidyverse)

source("Rlib/annotation.R")
cq<-paste0
ann=getAnnoteTable("Rlib/hgu133a.sqlite") %>%
    rownames_to_column("Probe") %>%
    select(Probe,gene_name,symbol) %>%
    distinct %>%
    tibble

tbl2=tibble(tbl) %>%
    left_join(eventTbl) %>%
    select(-Event) %>%
    rename(Event=Region,Gene=GENE) %>%
    left_join(geneAnn,by=c(Gene="SYMBOL")) %>%
    left_join(ann) %>%
    mutate(Description=ifelse(is.na(gene_name),GENENAME,gene_name)) %>%
    select(Band,Event,Gene,Description,Probe,Pvalue,log2FC,FDR,FC)

cghGeneTable=map(sheets,~readxl::read_xlsx(sTableFile,sheet=.)) %>% bind_rows

# Rewrite legacy HGNC aliases in cghGeneTable to current symbols so the join
# matches tbl2's u133a-derived symbols. Without this, e.g. SC5DL in
# cghGeneTable never matches SC5D in tbl2 and the row gets NA stats.
alias_map=AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys=unique(cghGeneTable$Gene),
        columns=c("SYMBOL"),
        keytype="ALIAS") %>%
    filter(!is.na(SYMBOL), ALIAS != SYMBOL) %>%
    distinct(ALIAS, .keep_all=TRUE) %>%
    rename(legacy=ALIAS, modern=SYMBOL)

cghGeneTable=cghGeneTable %>%
    left_join(alias_map, by=c(Gene="legacy")) %>%
    mutate(Gene=coalesce(modern, Gene)) %>%
    select(-modern)

tbl3=left_join(cghGeneTable,tbl2) %>%
    left_join(geneAnn,by=c(Gene="SYMBOL")) %>%
    mutate(Description=ifelse(is.na(Description),GENENAME,Description)) %>%
    select(-GENENAME) %>%
    arrange(Gene)

genes=tbl3 %>% filter(is.na(Description)) %>% distinct(Gene) %>% pull(Gene)
genes=intersect(genes, AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype="ALIAS"))

if(length(genes) > 0) {
    ann2=AnnotationDbi::select(
            org.Hs.eg.db::org.Hs.eg.db,
            keys=genes,
            columns=c("GENENAME"),keytype="ALIAS") %>%
        rename(SYMBOL=ALIAS) %>%
        distinct(SYMBOL,.keep_all=T) %>%
        arrange(SYMBOL) %>%
        tibble

    tbl3=tbl3 %>%
        left_join(ann2,by=c(Gene="SYMBOL")) %>%
        mutate(Description=ifelse(is.na(Description),GENENAME,Description)) %>%
        select(-GENENAME)
}

tbl3=tbl3 %>% group_split(Band)

names(tbl3)=map_vec(tbl3,~.$Band[1])
tbl3=tbl3 %>% map(~select(.,-Band)) %>% map(data.frame,check.names=F)
openxlsx::write.xlsx(tbl3,paste0("suppTable_",TABLE_NO,"_with_U133A__",SET_TAG,"_v3.xlsx"))

