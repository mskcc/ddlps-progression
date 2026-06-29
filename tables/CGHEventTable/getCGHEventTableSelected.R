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
    select(BIO_ID,everything())

openxlsx::write.xlsx(eventTblSel,"cghEventTable___Manuscript__SelectRegions1.xlsx")

# library(xlsx)
# write.csv(eventTable,file="cghEventTablePerBlockReducdedSet_v2.csv")

