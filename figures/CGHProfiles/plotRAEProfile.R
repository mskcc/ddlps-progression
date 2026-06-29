data(hg18)
data(cghUBM)
data(sampleTable)

##############
cArgs=commandArgs(trailing=T)
args=list(STYLE=1,TYPE="DD")
parseArgs=str_match(cArgs,"(.*)=(.*)")
apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
##
STYLE=as.numeric(args$STYLE)

wd.bioid=sampleTable$BIO_ID[sampleTable$aCGH=="Y" & sampleTable$TYPE=="WD"]
dd.bioid=sampleTable$BIO_ID[sampleTable$aCGH=="Y" & sampleTable$TYPE=="DD"]

breaks=ubm$breaks
offset=hg18$gOFFSET

breaks$Begin=breaks$Start+offset[breaks$Chr]
breaks$Stop=breaks$End+offset[breaks$Chr]

THRESHOLD=0.5
TYPE=args$TYPE

if(TYPE=="WD") {
    PLOTFILE=cc("raePlotA0D0_WD",paste0("s",STYLE),"v2")
    a0=ubm$a0[,wd.bioid]>THRESHOLD
    d0=ubm$d0[,wd.bioid]>THRESHOLD
} else if(TYPE=="DD") {
    PLOTFILE=cc("raePlotA0D0_DD",paste0("s",STYLE),"v2")
    a0=ubm$a0[,dd.bioid]>THRESHOLD
    d0=ubm$d0[,dd.bioid]>THRESHOLD
} else {
    cat("Invalid type",TYPE,"\n")
    quit()
}

xpos=(breaks$Begin+breaks$Stop)/2
meanA0=apply(a0,1,mean,na.rm=T)
meanD0=apply(-d0,1,mean,na.rm=T)

meanA0p=imputePolymorphicSomatic(ubm$poly,breaks,meanA0)
meanD0p=imputePolymorphicSomatic(ubm$poly,breaks,meanD0)

meanA0=meanA0p
meanD0=meanD0p

ampCol="#880000"
delCol="#000088"


if(STYLE==1 || STYLE==3) {
    png(file=cc(PLOTFILE,".png"),type="cairo",width=1200,height=900)
} else {
    pdf(file=cc(PLOTFILE,".pdf"),width=11,height=8)
}

plot(xpos,meanA0,type='n',ylim=c(-1,1),xlim=c(min(breaks$Begin),offset[23]),axes=F,xlab="",ylab="")
abline(h=0,lty=1,col=1)
abline(v=offset,lty=2,lwd=1,col=1)

if(STYLE==1) {
    for(i in seq(nrow(breaks)))
        rect((breaks$Begin[i]),0,(breaks$Stop[i]),meanA0[i],
            col=ampCol,border=ampCol)

    for(i in seq(nrow(breaks)))
        rect((breaks$Begin[i]),0,(breaks$Stop[i]),meanD0[i],
            col=delCol,border=delCol)

} else if (STYLE==2) {

    ii=seq(nrow(breaks))[-1]
    for(i in ii) {
        xx=c(breaks$Stop[i-1],breaks$Begin[i],breaks$Stop[i])
        yy=c(meanA0[i-1],meanA0[i],meanA0[i])
        lines(xx,yy,col=ampCol,lwd=3)
    }

    for(i in ii) {
        xx=c(breaks$Stop[i-1],breaks$Begin[i],breaks$Stop[i])
        yy=c(meanD0[i-1],meanD0[i],meanD0[i])
        lines(xx,yy,col=delCol,lwd=3)
    }
} else if (STYLE==3) {
    points(xpos,meanA0p,type='s',col=ampCol,lwd=2)
    points(xpos,meanD0p,type='s',col=delCol,lwd=2)
}

chrL=(offset[-len(offset)]+diff(offset)/2)[1:22]
text(chrL,-.75+.05*(seq(chrL)%%2),seq(chrL))


box()
axis(2)

dev.off()

