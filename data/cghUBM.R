#
# Limit to reduced set of samples
#

source("R/helpers.R")
source("tools.R")
reducedSet=scan("db/cghCBSRdataFiles","")
reducedSet=fixCGHSlideNames(reducedSet)

source("crdb.R")
crdb=crdb[crdb$ASSAY=="aCGH",]
crdb$CTAGS=fixCGHSlideNames(crdb$FILE_NAME)

load("db/CragoProgressionC-lesions.Rdata")
ubm=list()
ubm$breaks=data.frame(lesions$breaks)

raeTags=fixCGHSlideNames(colnames(lesions$a0))
reducedII=(raeTags %in% reducedSet)

ubm$a0=lesions$a0[,reducedII]
colnames(ubm$a0)=crdb$BIO_ID[match(raeTags[reducedII],crdb$CTAGS)]

ubm$d0=lesions$d0[,reducedII]
colnames(ubm$d0)=crdb$BIO_ID[match(raeTags[reducedII],crdb$CTAGS)]

ubm$s0=ifelse(ubm$d0>ubm$a0,-ubm$d0,ubm$a0)

#
# Check if we have the poly sites for this
# UBM model
#

require(digest)
signature=digest(c(ubm$breaks$Chr,ubm$breaks$Start,ubm$breaks$End))
polyFile=cc("db/poly",signature,".Rdata")

if(!file.exists(polyFile)){
    cat("Computing polymorphic sites, will take awhile...")
    poly=getPolymorphicCNV(ubm$breaks)
    save(poly,file=polyFile,compress=T)
    cat("\n")
} else {
    load(polyFile)
}
ubm$poly=poly

rm("poly")
rm("crdb")
rm("CRDB.FILE")
rm("crdb.full")
rm("crdb.orig")
rm("lesions")
rm("lipoTypes")
rm("raeTags")
rm("reducedII")
rm("reducedSet")

