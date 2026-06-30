source("R/helpers.R")

data(hg18)
data(cghUBM)
data(sampleTable)


wd.bioid=sampleTable$BIO_ID[sampleTable$aCGH=="Y" & sampleTable$TYPE=="WD"]
dd.bioid=sampleTable$BIO_ID[sampleTable$aCGH=="Y" & sampleTable$TYPE=="DD"]

breaks=ubm$breaks
offset=hg18$gOFFSET

breaks$Begin=breaks$Start+offset[breaks$Chr]
breaks$Stop=breaks$End+offset[breaks$Chr]

data(regions)
rii=3
data(u133a)
#genes=c("TCF21","SASH1","UST","TAB2")
genes=c("SASH1","CCDC28A","TAB2")
genes2=scan("geneSet2Chr6","")

xpos=(breaks$Begin+breaks$Stop)/2

THRESHOLD=0.5
a0=ubm$a0[,wd.bioid]>THRESHOLD
meanA0WD=apply(a0,1,mean,na.rm=T)
meanA0WD=imputePolymorphicSomatic(ubm$poly,breaks,meanA0WD)

a0=ubm$a0[,dd.bioid]>THRESHOLD
meanA0DD=apply(a0,1,mean,na.rm=T)
meanA0DD=imputePolymorphicSomatic(ubm$poly,breaks,meanA0DD)

ii=which(regionsM$CHR[rii]==breaks$Chr & regionsM$START[rii]<=breaks$Start & regionsM$END[rii]>=breaks$End)

ampColDD=rgb(1,5/16,5/16)
ampColWD=rgb(10/16,3/16,3/16)

###################################################################################
labelGenes<-function(genes){
    for(i in seq(len(genes))){
        gg=genes[i]
        print(gg)
        probe=rownames(u133a$annote)[u133a$annote$SYMBOL==gg][1]
        ggS=u133a$annoteHG18[probe,"START"]
        ggE=u133a$annoteHG18[probe,"END"]
        rect(ggS,0,ggE,1,col=8,border=NA)
    }

    cat("USR=",paste(par()$usr),"\n\n")

    dLabel=-240000
    for(i in seq(len(genes))){
        gg=genes[i]
        print(gg)
        probe=rownames(u133a$annote)[u133a$annote$SYMBOL==gg][1]
        ggS=u133a$annoteHG18[probe,"START"]
        ggE=u133a$annoteHG18[probe,"END"]
        abline(v=(ggS+ggE)/2,col=1,lwd=1,lty=2)
        text((ggS+ggE)/2+dLabel,1.5*par()$usr[3],gg,pos=1,srt=45,xpd=T,cex=0.8)
    }
}

###################################################################################
pdf(file="chr6RegionTest_v3.pdf",width=11,height=8.5)

plot(xpos,meanA0DD,type='n',ylim=c(0,.35),
    xlim=c(regionsM$START[rii], regionsM$END[rii]),
    axes=F,xlab="",ylab="")

abline(h=0,lty=2,col=8)

labelGenes(genes)

for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0DD[i-1],meanA0DD[i],meanA0DD[i])
    lines(xx,yy,col=2,lwd=3)
}


for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0WD[i-1],meanA0WD[i],meanA0WD[i])
    lines(xx,yy,col=ampColWD,lwd=3)
}


box()
axis(2)
legend(par()$usr[2]*.9,.325,c("DD","WD"),lwd=2,col=c(ampColDD,ampColWD))
dev.off()

###################################################################################
###################################################################################
pdf(file="chr6RegionTest_WDOnly_v3.pdf",width=11,height=8.5)

plot(xpos,meanA0WD,type='n',ylim=c(0,.2),
    xlim=c(regionsM$START[rii], regionsM$END[rii]),
    axes=F,xlab="",ylab="",main="Fraction Amp [WDonly]")

abline(h=0,lty=2,col=8)

labelGenes(genes)

for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0WD[i-1],meanA0WD[i],meanA0WD[i])
    lines(xx,yy,col=ampColWD,lwd=3)
}

box()
axis(2)

plot(xpos,meanA0WD,type='n',ylim=c(0,.2),
    xlim=c(regionsM$START[rii], regionsM$END[rii]),
    axes=F,xlab="",ylab="",main="Fraction Amp [WDonly]")

abline(h=0,lty=2,col=8)

labelGenes(genes2)

for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0WD[i-1],meanA0WD[i],meanA0WD[i])
    lines(xx,yy,col=ampColWD,lwd=3)
}

box()
axis(2)

plot(xpos,meanA0WD,type='n',ylim=c(0,.2),
    xlim=c(regionsM$START[rii], regionsM$END[rii]),
    axes=F,xlab="",ylab="",main="Fraction Amp [WDonly]")

abline(h=0,lty=2,col=8)

for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0WD[i-1],meanA0WD[i],meanA0WD[i])
    lines(xx,yy,col=ampColWD,lwd=3)
}

box()
axis(2)

dev.off()

##########################################################################
##########################################################################

pdf(file="chr6Region_WD+DD_v4.pdf",width=11,height=8.5)

plot(xpos,meanA0DD,type='n',ylim=c(0,.35),
    xlim=c(regionsM$START[rii], regionsM$END[rii]),
    axes=F,xlab="",ylab="")

abline(h=0,lty=2,col=8)

#labelGenes(genes)

for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0DD[i-1],meanA0DD[i],meanA0DD[i])
    lines(xx,yy,col=2,lwd=3)
}


for(i in ii) {
    xx=c(breaks$End[i-1],breaks$Start[i],breaks$End[i])
    yy=c(meanA0WD[i-1],meanA0WD[i],meanA0WD[i])
    lines(xx,yy,col=4,lwd=3)
}


box()
axis(2)
xLegend=par()$usr[2]-.125*diff(par()$usr[1:2])
legend(xLegend,.335,c("DD","WD"),lwd=3,col=c(2,4))
dev.off()

