regions=scan("db/survialRegions","")


regionsM=data.frame(
    EVENT=sapply(strsplit(regions,"\\."),"[",1),
    CHR=sapply(strsplit(regions,"\\."),"[",2),
    START=as.numeric(sapply(strsplit(regions,"\\."),"[",3)),
    END=as.numeric(sapply(strsplit(regions,"\\."),"[",4))
    )

