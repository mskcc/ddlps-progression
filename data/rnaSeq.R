rnaSeq=list()
dx=read.delim("db/Proj_3704_Merge_GeneCounts.txt.gz")
i.pc=dx$CLASS=="protein_coding"
ds=dx[i.pc,3:ncol(dx)]
rownames(ds)=dx[i.pc,1]

rnaSeq$ds=ds
rnaSeq$dFull=dx

rm(ds)
rm(dx)
rm(i.pc)
