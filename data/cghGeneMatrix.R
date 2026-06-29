source("tools.R")

geneCGH=read.delim("db/GeneMatrix.txt",row.names=1)
colnames(geneCGH)=fixCGHSlideNames(colnames(geneCGH))

reducedSet=scan("db/cghCBSRdataFiles","")
reducedSet=fixCGHSlideNames(reducedSet)

source("crdb.R")
crdb=crdb[crdb$ASSAY=="aCGH",]
crdb$CTAGS=fixCGHSlideNames(crdb$FILE_NAME)

reducedSetII=(colnames(geneCGH) %in% reducedSet)

geneCGH=geneCGH[,reducedSetII]
colnames(geneCGH)=crdb$BIO_ID[match(colnames(geneCGH),crdb$CTAGS)]

rm("crdb")
rm("CRDB.FILE")
rm("crdb.full")
rm("crdb.orig")
rm("fixCGHSlideNames")
rm("getPolymorphicCNV")
rm("imputePolymorphicSomatic")
rm("lipoTypes")
rm("raeCNV")
rm("reducedSet")
rm("reducedSetII")
