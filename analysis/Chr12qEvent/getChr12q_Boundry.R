data(cghUBM)
data(progressionSet)

ii.12=ubm$breaks$Chr==12

# Threshold A0 to .25 and then find regions where mean
# over the 8 Progression First WD >0.95 (ie all of them)
#

firstA0=apply(ubm$a0[ii.12,progressionSet$FIRST]>.25,1,mean)
segments=ubm$breaks[ii.12,][ which(firstA0>.95),]


mergedSegs=NULL
merge=NULL
for(i in seq(nrow(segments)-1)) {
    if(segments$StartIdx[i+1]==segments$EndIdx[i]+1){
        segments$Start[i+1]=segments$Start[i]
        segments$Size[i+1]=segments$Size[i+1]+segments$Size[i]
        segments$Length[i+1]=segments$Length[i+1]+segments$Length[i]
        segments$StartIdx[i+1]=segments$StartIdx[i]
        segments[i,]=rep(-1,ncol(segments))
    }
}

event=segments[segments$Chr!=-1,]
write.xls(event,"chr12qEvent.txt")
write.csv2(event,"data/chr12qEvent.csv",row.names=F)
