source("funcs.R")
data(sampleTable)
data(cghEventsPerBlock)

require(tidyverse)
require(bedr)
require(openxlsx)

aCGH.SampleTable=tibble(sampleTable[sampleTable$aCGH=="Y",])

eventTableLong=cghEvents %>%
    rownames_to_column("BIO_ID") %>%
    tibble %>%
    filter(BIO_ID %in% aCGH.SampleTable$BIO_ID) %>%
    gather(Event,CNV,-BIO_ID) %>%
    left_join(aCGH.SampleTable)


evtTbl=map(unique(eventTableLong$Event),parseEventTag) %>%
    map(as_tibble) %>%
    bind_rows %>%
    mutate.getRegion

bands=read_csv("table4Bands.csv") %>% mutate.getRegion %>% mutate(TTag=paste(BAND,Call,WD,DD,sep=":"))

evt.regions=bedr.sort.region(unique(evtTbl$Region))
band.regions=bedr.sort.region(unique(bands$Region))

evt2bands=bedr(
    input=list(b=band.regions,a=evt.regions),
    method="intersect",
    params="-loj -sorted") %>%
    tibble %>%
    filter(V4!=".") %>%
    mutate(BandRegions=getRegion(V4,V5,V6)) %>%
    rename(EvtRegions=index) %>%
    select(matches("Regions")) %>%
    left_join(bands,by=c(BandRegions="Region")) %>%
    select(Region=EvtRegions,BAND,Chr,Start,Stop) %>%
    right_join(evtTbl,by="Region") %>%
    select(Event,BAND)

ee=eventTableLong %>%
    left_join(evt2bands) %>%
    filter(!is.na(BAND)) %>%
    count(BAND,Event,CNV,TYPE) %>%
    filter(CNV=="X") %>%
    spread(TYPE,n) %>%
    mutate(Call=gsub("\\..*","",Event)) %>%
    mutate(TTag=paste(BAND,Call,WD,DD,sep=":")) %>%
    left_join(bands) %>%
    filter(!is.na(Region))

featTbl=read_tsv("data/raw/CGH/rae/FEAT.file",comment="#") %>%
    filter(Type!="Cytoband") %>%
    rename(Stop=End) %>%
    mutate(Chr=paste0("chr",Chr)) %>%
    mutate.getRegion %>%
    select(Name,Region)

ee.regions=bedr.sort.region(ee$Region)
feat.regions=bedr.sort.region(unique(featTbl$Region))

feat2bands =
    bedr(input=list(a=ee.regions,b=feat.regions),method="intersect",params="-loj -sorted -F 0.1") %>%
    tibble %>%
    mutate(Feat.Region=getRegion(V4,V5,V6)) %>%
    rename(Band.Region=index) %>%
    select(matches("Region")) %>%
    left_join(featTbl,by=c(Feat.Region="Region")) %>%
    left_join(ee,by=c(Band.Region="Region")) %>%
    select(BAND,Event,Name,DD,WD) %>%
    left_join(featTbl,by="Name") %>%
    mutate(Pos=as.numeric(str_match(Region,":(.*)-")[,2])) %>%
    arrange(BAND,Pos)

ll=feat2bands %>% mutate(Gene=gsub(".*:","",Name)) %>% select(BAND,Event,Gene,WD,DD) %>% distinct %>% group_split(BAND)
names(ll)=map(ll,"BAND") %>% map(unique) %>% unlist
write.xlsx(ll,"raeGeneTableSigRegionsV1.xlsx")

#ite.xlsx(ll,)
