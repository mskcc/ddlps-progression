#library(siggenes)

P.filter <- function(calls,fact,cutP) {
  numP = tapply(calls=="P",fact,sum)
  return(any(numP>cutP))
}

filter.P <- function(dp,cl) {
  cutP = table(cl)/4
  ig=apply(dp,1,P.filter,cl,cutP)
  return(ig)
}

mk.diff.table <- function(ds,cl,gg,pv,sampA=F,sampB=F) {

  mm=t(apply(ds[gg,],1,function(x,ff){tapply(x,ff,mean)},cl))
  fc=apply(mm,1,function(x){log2FC(x[1]-x[2])})
  
  ans=data.frame(p.value=pv,2^mm,fc)
  if(is.logical(sampA) && !sampA) {
    sampA = levels(factor(cl))[1]
    sampB = levels(factor(cl))[2]
  }
  
  mm.lab=c(paste("mean",sampA,sep="."),paste("mean",sampB,sep="."))
  fc.lab=paste("FC",sampA,sampB,sep=".")
  colnames(ans)[2:4]=c(mm.lab,fc.lab)

  return(ans)
}

##
## Fix for siggenes 1.2.17 (BioC 1.6)
##
sam.diff.table <- function(ds,cl,res,sampA="A",sampB="B") {

  gg=rownames(res@mat.sig)
  pv=res@mat.sig$p.value

  ans=mk.diff.table(ds,cl,gg,pv,sampA,sampB)

  return(ans)

}

calc.bayes.qv <- function(tt) {

  tstat=tt$t.score
  df=tt$dof1
  pval <- 2 * pt(-abs(tstat), df)
  qv=qvalue.cal(pval,1)

  names(qv)=rownames(tt)
  names(pval)=rownames(tt)

  qv=sort(qv)
  return(list(pv=pval,qv=qv))

}

log2FC <- function(x) {
  return(ifelse(x<0,-2^(-x),2^x))
}

rn <- function(x, digits = 3){
  if (is.null(x)) 
    NULL
  else {
    if (is.matrix(x) && ncol(x) == 1) 
      x <- x[, 1]
    round(x, digits = digits)
  }
}



glog2 <- function(x,p0=0,p1=1) {
  return((asinh(p0+p1*x)-log(2*p1))/log(2))
}

#
# Only for p0=0, p1=1
#
gexp2 <- function(x){
  return(sinh(log(2)*x+log(2)))
}

#
# abs Derivative
#
adgexp2 <- function(x) {
  return(abs(log(2)*cosh(log(2)*x+log(2))))
}
