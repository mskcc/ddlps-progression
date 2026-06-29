parseEventTag<-function(ee) {

    pp=str_match(ee,"(Amp|Del)\\.(\\d+)\\.(\\d+)\\.(\\d+)")
    names(pp)=c("Event","Call","Chr","Start","Stop")
    pp[3]=paste0("chr",pp[3])
    as.list(pp)
}

getRegion<-function(chr,start,stop) {
    paste0(chr,":",start,"-",stop)
}

mutate.getRegion<-function(.data,...) {
    mutate(.data,Region=getRegion(Chr,Start,Stop))
}
