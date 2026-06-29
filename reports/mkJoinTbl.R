require(readxl)
require(tidyverse)
tbl2=read_xlsx("sampleTableProgressionV10.xlsx",sheet=2)
clin=read_xlsx("sampleTableProgressionV10.xlsx",sheet=3)

clinR=clin %>%
    select(P_MRN,BIO_ID,matches("(^SARCOMA|^MOLD|DATE|^PART)"),BOX,FREEZER,RACK) %>%
    distinct
jTbl=left_join(tbl2,clinR,by="BIO_ID")
require(openxlsx)
write.xlsx(as.data.frame(jTbl),file="sampleTableProgressionV10__ClinJoin.xlsx",rowNames=FALSE,keepNA=FALSE)
