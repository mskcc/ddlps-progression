data(sampleTable)
data(cghEventsPerBlock)

aCGH.SampleTable=sampleTable[sampleTable$aCGH=="Y",]
eventTable=data.frame(
    aCGH.SampleTable[,c(2,3,4,5,7)],
    cghEvents[aCGH.SampleTable$BIO_ID,])

require(tidyverse)
regions=scan("manuscriptCGHRegions_2023-02-25.txt","")
eventTblSel=tibble(eventTable) %>%
    select(1:5,all_of(regions)) %>%
    left_join(tibble(sampleTable) %>% filter(aCGH=="Y") %>%
    select(BIO_ID,P_MRN,MOLD_DATE,PROCUREMENT_DATE,TYPE)) %>%
    select(BIO_ID,everything()) %>%
    gather(Event,Call,matches("^(Amp|Del)\\.")) %>%
    filter(Call=="X") %>%
    separate(Event,c("EType","Chr","Start","End"),remove=F) %>%
    mutate(Chr=paste0("chr",Chr)) %>%
    type_convert %>%
    mutate(UUID=row_number())

xref=read_tsv("ucsc_hg18__kgXref.txt.gz") %>%
    rename(X1=`#kgID`) %>%
    filter(!is.na(refseq) & !is.na(protAcc)) %>%
    select(X1,geneSymbol,description)

genes=read_tsv("ucsc_hg18__knownGene.txt.gz",col_names=F) %>%
    left_join(xref) %>%
    select(Chr=X2,Start=X4,End=X5,X1,geneSymbol,description) %>%
    type_convert

require(tidygenomics)

valid_gene_overlaps=eventTblSel %>%
    select(UUID,Chr,Start,End) %>%
    genome_intersect(select(genes,X1,Chr,Start,End)) %>%
    left_join(select(eventTblSel,UUID,E.Start=Start,E.End=End,EType),by="UUID") %>%
    left_join(select(genes,X1,G.Start=Start,G.End=End),by="X1") %>%
    mutate(LEN=End-Start,GLEN=G.End-G.Start,Pog=LEN/GLEN) %>%
    filter(Pog==1)

tbl2=valid_gene_overlaps %>%
    select(UUID,X1) %>%
    left_join(eventTblSel) %>%
    left_join(genes,by="X1") %>%
    filter(!is.na(geneSymbol)) %>%
    distinct(UUID,geneSymbol,.keep_all=T) %>%
    select(geneSymbol,description,Event,TYPE,BIO_ID,Call) %>%
    group_split(TYPE) %>%
    map(\(x){
        distinct(x) %>% spread(BIO_ID,Call) %>% arrange(Event)
    })

names(tbl2)=map(tbl2,\(x) x$TYPE[1]) %>% unlist

openxlsx::write.xlsx(tbl2,"cghGeneTable___Manuscript___SelectRegions1.xlsx")
