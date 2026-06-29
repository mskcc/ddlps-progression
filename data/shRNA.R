dd.shRNA=read.delim("db/averageSignal_3Sort_M3_.txt")
colnames(dd.shRNA)[1]="SYMBOL"
rownames(dd.shRNA)=dd.shRNA$SYMBOL
dd.shRNA$RANK=seq(nrow(dd.shRNA))
dd.shRNA$PCT.RANK=100*seq(nrow(dd.shRNA))/nrow(dd.shRNA)
dat.cols=grep("_SOFT_TISSUE",colnames(dd.shRNA))
dd.shRNA$MIN.KD=apply(dd.shRNA[,dat.cols],1,min)
dd.shRNA$MED.KD=apply(dd.shRNA[,dat.cols],1,median)

shRNA=list()
shRNA$dd=dd.shRNA
shRNA$file="db/averageSignal_3Sort_M3_.txt"
rm(dd.shRNA)
