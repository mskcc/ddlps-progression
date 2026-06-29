##
## Sat May 17 22:15:25 EDT 2008
##

library(RSQLite)

dbDir="~/Work/AnnoteDBs"
dbFiles=paste(dbDir,dir(dbDir),sep="/")

getAnnoteTable <- function(dbFile) {
  drv=dbDriver("SQLite")
  dbList=dir(dbDir)
  con=dbConnect(drv,dbFile)
  probes=dbReadTable(con,"probes",row.names="probe_id")
  probeIDs=as.character(probes[,2])
  ##genes=dbReadTable(con,"genes",row.names="X_id")
  info=dbReadTable(con,"gene_info")
  rownames(info)=info$X_id
  ann=cbind(probes[,1,drop=F],info[probeIDs,])
  maploc=dbReadTable(con,"chromosome_locations")
  maploc=maploc[!duplicated(maploc[,"X_id"]),]
  rownames(maploc)=maploc$X_id
  maploc$Pos=cq("Chr",maploc[,2],":",abs(maploc[,3]))
  ann=cbind(ann,maploc[probeIDs,"Pos",drop=F])
  ii=which(is.na(ann[,1]))
  if(length(ii)>0)
    ann[ii,1]=cc("_id",rownames(ann)[ii])
  ii=which(is.na(ann[,2]))
  if(length(ii)>0) {
    ann[ii,"symbol"]=ann[ii,1]
    ann[ii,"gene_name"]=ann[ii,1]
  }

  dbDisconnect(con)
  dbUnloadDriver(drv)
  return(ann)
}

load.GO.ann <- function(dbFile) {
  drv=dbDriver("SQLite")
  dbList=dir(dbDir)
  con=dbConnect(drv,dbFile)
  go=list()
  go$bp=dbReadTable(con,"go_bp_all")
  go$mf=dbReadTable(con,"go_mf_all")
  go$cc=dbReadTable(con,"go_cc_all")
  cgo=dbConnect(drv,"~/Work/AnnoteDBs/GO.sqlite")
  go$terms=dbReadTable(cgo,"go_term",row.names="go_id")
  return(go)
}

get.go <- function(gene.ids,go)
{
  ans=list()
  ans$X_id=sapply(gene.ids,function(x){unique(ann[rownames(ann)==x,"X_id"])})
  ans$GO.BP=sapply(ans$X_id,function(x){
    paste(unique(go$terms[go$bp[go$bp$X_id==x,"go_id"],"term"]),collapse=" // ")})
  ans$GO.MF=sapply(ans$X_id,function(x){
    paste(unique(go$terms[go$mf[go$mf$X_id==x,"go_id"],"term"]),collapse=" // ")})
  ans$GO.CC=sapply(ans$X_id,function(x){
    paste(unique(go$terms[go$cc[go$cc$X_id==x,"go_id"],"term"]),collapse=" // ")})
  return(data.frame(ans))
}

##
##
##
