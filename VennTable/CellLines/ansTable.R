ansTable<-function (fit, coef = 1, q.cut=0.05, lfc=0, number = -1,
                      ann=NULL,vstParam=NULL,sig.dig=2)
{
  if(number<0) number=nrow(fit$p.value)
  logFC <- as.matrix(fit$coefficients)[, coef]
  FC <- log2FC(logFC)
  p.value <- as.matrix(fit$p.value)[, coef]
  B <- as.matrix(fit$lods)[, coef]
  FDR <- p.adjust(p.value, method = 'fdr')
  ord=order(p.value,decreasing=FALSE)
  if (q.cut < 1 | lfc > 0) {
    sig <- (FDR < q.cut) & (abs(logFC) >= lfc)
    if (any(is.na(sig)))
      sig[is.na(sig)] <- FALSE
    nsig <- sum(sig)
    if (nsig == 0)
      return(data.frame())
    top <- ord[sig[ord]]
    if (number < nsig)
      top <- top[1:number]
  }
  else {
    top <- ord[1:number]
  }
  tab <- data.frame(FDR=(FDR[top]),
                    FC = round(FC[top],sig.dig))

  if (!is.null(fit$ds)) {
    levels=colnames(fit$design)
    type=character(nrow(fit$design))
    for(ii in seq(ncol(fit$design)))
      type[which(fit$design[,ii]!=0)]=levels[ii]
    type=factor(type)
    avgExp=t(apply(fit$ds[top,,drop=F],1,function(x){tapply(x,type,mean)}))
    if(!is.null(vstParam)) {
      avgExp=inverseVST(avgExp,parameter=vstParam)
    } else {
      avgExp=2^avgExp
    }
    typesUsed=names(which(fit$contrast[,coef]!=0))
    tab <- data.frame(tab, avg=round(avgExp[,typesUsed,drop=F],sig.dig))
  }

  if(!is.null(ann)) {
    pID=as.character(rownames(tab))
    tab=data.frame(
      SYM=ann[pID,"symbol"],
      tab,
      NAME=ann[pID,"gene_name"])
  }

  return(tab)
}






