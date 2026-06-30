source("R/helpers.R")
require(limma)
source("diffGenes.R")

data(cellLines)

ds=cellLines$ds
key=cellLines$key
type=factor(key$Class)

design=model.matrix(~0+type)
colnames(design)=levels(type)
fit=lmFit(ds,design)

col.design=colnames(design)
N.t=length(col.design)
cm=matrix(0,nrow=N.t,ncol=N.t*(N.t-1))
dimnames(cm)=list(
        Levels=colnames(design),
        Contrasts=seq(ncol(cm)))

#######################################
#### Customize this to select the
#### correct (wanted) contrasts
#######################################

select.cont <- function(c.A,c.B) {
  if(c.B %in% c("Control","Untreated")) {
    if(! (c.A %in% c("Control","Untreated"))) {
    return(TRUE)
    }
  }
  return(FALSE)
}

#######################################
##
## All pairwise contrasts
##

con=1
for(i in 1:N.t) {
  for(j in 1:N.t) {
    if(select.cont(col.design[i],col.design[j])) {
      cm[col.design[i],con]=1
      cm[col.design[j],con]=-1
      colnames(cm)[con]=paste(col.design[i],col.design[j],sep="-")
      con = con + 1
    }
  }
}

cm.o=cm[,1:(con-1),drop=FALSE]
cm=cm.o
cnames=colnames(cm.o)

#
# For paper only want
#    CDK4_inhib-Untreated
#    MDM2_KD_380-Control
#

cm=cm[,c("CDK4_inhib-Untreated","MDM2_KD_380-Control"),drop=F]

fit2=contrasts.fit(fit,cm)
fit2=eBayes(fit2)

source("ansTable.R")
source("fileTools.R")
require(lumiHumanAll.db)
require(annotate)

ALL=TRUE
if(ALL==TRUE) {
  q.cut=1.01
  lfc=log2(1)
} else {
  q.cut=0.05
  lfc=log2(2)
}

fit2$ds=ds

ans=list()
for(ci in seq(ncol(fit2$contrasts))) {
  print(ci)
  try({ans[[colnames(fit2$contrasts)[ci]]]=ansTable(fit2,coef=ci,q.cut=q.cut,lfc=lfc)})
}

for(ii in seq(ans)) {
  if(nrow(ans[[ii]])>0) {
    gid=rownames(ans[[ii]])
    aa=data.frame(
      SYMBOL=unlist(lookUp(gid,"lumiHumanAll.db","SYMBOL")),
      ans[[ii]],
      NAME=unlist(lookUp(gid,"lumiHumanAll.db","GENENAME")))

    file=cc("diffGenes",DATE(),names(ans)[[ii]],'FDR',q.cut,"FC",2^lfc,".txt")
    print(file)
    write.xls(aa,file)
  }
}



