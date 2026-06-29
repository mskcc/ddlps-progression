data(u133a)
data(sampleTable)
ii.u133.nfwd=sampleTable$U133A=="Y" & sampleTable$TYPE %in% c("NF","WD")
samples=sampleTable$BIO_ID[ii.u133.nfwd]

ds=u133a$ds[,samples]
class=sampleTable$TYPE[ii.u133.nfwd]
bioids=sampleTable$BIO_ID[ii.u133.nfwd]
type=factor(class)

joinTbl=read.delim("joinTableCragoProgression_v14_.txt")
probes=joinTbl$PROBE.U133A

dds=ds[probes,]

dZ=sweep(dds,1,apply(dds,1,mean),"-")
dZZ=sweep(dZ,1,apply(dZ,1,sd),"/")

library(gplots)


dZZc=dZZ
dZZc[dZZc>4]=4
dZZc[dZZc< -4]=-4

classCols=c("#0000EE","#EE0000")
rc=classCols[factor(class)]

pdf(file="heatMap_v14.pdf",width=8.5,height=11)

heatmap.2(t(dZZc),
    trace="none",
    dendro="row",
    col=colorpanel(32,"blue","white","red"),
    symbreaks=TRUE,labRow=bioids,
    labCol="",
    RowSideColors=rc)

heatmap.2(t(dZZc),
    trace="none",
    dendro="row",
    col=colorpanel(32,"blue","white","red"),
    symbreaks=TRUE,labRow=class,
    labCol="",
    RowSideColors=rc)

dev.off()

#
# === CLASSIFICATION ACCURACY OF THE HEATMAP (Figure 1A) ===
# The published Figure 1A reads three brackets off the row dendrogram
# ("21/22 normal fat", "35/42 WDLS", "26/26 WDLS"). Reproduce that here so the
# reported numbers come from the pipeline, not from eyeballing the figure.
#
# The dendrogram is exactly the one heatmap.2 drew above: it clusters the rows
# of t(dZZc) (= the samples) with heatmap.2's defaults, distfun=dist (Euclidean)
# and hclustfun=hclust (complete linkage). We re-run those defaults explicitly,
# cut at k=3 (the three brackets in the figure), label each bracket by its
# majority TYPE, and compute accuracy / sensitivity / specificity / precision
# treating WDLS (tumor) as the positive class.
#
# See methods_heatmap.md for the full rationale (terminology + metric choice).
#
hc=hclust(dist(t(dZZc)))            # heatmap.2 defaults: Euclidean + complete
grp=cutree(hc,k=3)                  # k=3 = the three dendrogram brackets

tab3=table(bracket=grp,truth=class)
maj=apply(tab3,1,function(r) colnames(tab3)[which.max(r)])   # majority TYPE per bracket
pred=factor(maj[as.character(grp)],levels=c("WD","NF"))
truth=factor(class,levels=c("WD","NF"))

TP=sum(pred=="WD" & truth=="WD")    # positive = WDLS (tumor)
FP=sum(pred=="WD" & truth=="NF")
FN=sum(pred=="NF" & truth=="WD")
TN=sum(pred=="NF" & truth=="NF")
N=TP+FP+FN+TN

accuracy   =(TP+TN)/N
sensitivity=TP/(TP+FN)              # recall / TPR; WDLS detected
specificity=TN/(TN+FP)             # TNR; normal fat detected
precision  =TP/(TP+FP)             # PPV

con=file("heatmap_accuracy_v14.txt","w")
writeLines(c(
  "Figure 1A heatmap -- classification accuracy",
  "Clustering: Euclidean distance + complete linkage (heatmap.2 defaults), cut at k=3.",
  "Positive class = WDLS (tumor); Negative = normal fat (NF).",
  "",
  sprintf("Entities clustered: %d  (truth: %d WD, %d NF)",N,sum(truth=="WD"),sum(truth=="NF")),
  "",
  "Three dendrogram brackets (top -> bottom in figure):"),con)
for(cl in unique(grp[hc$order])){
    inb=grp==cl
    nNF=sum(inb & class=="NF"); nWD=sum(inb & class=="WD")
    m=ifelse(nWD>=nNF,"WD","NF"); corr=ifelse(m=="WD",nWD,nNF)
    writeLines(sprintf("  bracket n=%2d: NF=%2d WD=%2d -> called %s (%d/%d correct)",
        sum(inb),nNF,nWD,m,corr,sum(inb)),con)
}
writeLines(c(
  "",
  "Confusion matrix (positive = WD):",
  sprintf("  TP=%d  FP=%d  FN=%d  TN=%d",TP,FP,FN,TN),
  "",
  sprintf("Accuracy    = (TP+TN)/N   = %d/%d = %.1f%%",TP+TN,N,100*accuracy),
  sprintf("Sensitivity = TP/(TP+FN)  = %d/%d = %.1f%%   (WDLS recall)",TP,TP+FN,100*sensitivity),
  sprintf("Specificity = TN/(TN+FP)  = %d/%d = %.1f%%   (normal-fat recall)",TN,TN+FP,100*specificity),
  sprintf("Precision   = TP/(TP+FP)  = %d/%d = %.1f%%   (PPV for WDLS)",TP,TP+FP,100*precision)),con)
close(con)

cat(readLines("heatmap_accuracy_v14.txt"),sep="\n")
cat("\n")


