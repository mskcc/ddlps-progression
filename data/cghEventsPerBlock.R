source("tools.R")

cghEvents=read.delim("db/CragoProgressionC_EventMatrixV2_20160419_.txt",row.names=1)
cghTags=fixCGHSlideNames(rownames(cghEvents))

#
# Limit to reduced set of samples
#

reducedSet=scan("db/cghCBSRdataFiles","")
reducedSet=fixCGHSlideNames(reducedSet)
ii=cghTags %in% reducedSet
cghTags=cghTags[ii]
cghEvents=cghEvents[ii,]

source("crdb.R")
crdb=crdb[crdb$ASSAY=="aCGH",]
crdb$CTAGS=fixCGHSlideNames(crdb$FILE_NAME)

rownames(cghEvents)=crdb$BIO_ID[match(cghTags,crdb$CTAGS)]

rm("cghTags")
rm("crdb")
rm("CRDB.FILE")
rm("crdb.full")
rm("crdb.orig")
rm("lipoTypes")
rm("reducedSet")
