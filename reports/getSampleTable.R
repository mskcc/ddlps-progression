TABLENAME="sampleTableProgressionV10"

require(openxlsx)
require(knitr)

data(crdb)

xx=xtabs(~ TYPE + fASSAY, data=crdb[crdb$REDUCED.SET=="Y",])
xx=as.data.frame.matrix(xx)
xx=rbind(xx,TOTAL=colSums(xx))

tbl1=xx

data(sampleTable)

mafSamples=sampleTable[sampleTable$CustomArray=="X",]

# Patient Level Counts
xx=table(sapply(strsplit(unique(paste(mafSamples$P_MRN,mafSamples$TYPE))," "),"[",2))
xx=xx[c("WD","DD")]
mafCounts=data.frame(Patient.Custom=as.numeric(xx))
rownames(mafCounts)=names(xx)
xy=table(mafSamples$TYPE)
xy=xy[c("WD","DD")]
mafCounts=cbind(mafCounts,data.frame(Sample.Custom=as.numeric(xy)))
mafCounts=rbind(mafCounts,colSums(mafCounts))
rownames(mafCounts)[3]="TOTAL"
mafCounts=rbind(c(0,0),mafCounts)
rownames(mafCounts)[1]="NF"
xx=cbind(tbl1,mafCounts)

write.xlsx(
    list(
        Summary    = as.data.frame.matrix(xx),
        AssayTable = sampleTable,
        Full.CRDB  = crdb,
        Excluded   = crdb.orig[crdb.orig$EXCLUDED=="Y",]
    ),
    paste0(TABLENAME,".xlsx"),
    rowNames=c(TRUE, FALSE, FALSE, FALSE)
)

