hg18=read.delim("db/human.hg18.genome",header=F,col.names=c("CHROM","LENGTH"))
hg18$gOFFSET=cumsum(c(0,hg18$LENGTH[-nrow(hg18)]))
