CRDB.FILE="raw/CRDB/FullDataExport_Nick_Socci_04_05_2016___PaperFREEZE_2017_05_03.txt"
crdb.orig=read.delim(CRDB.FILE)
crdb=crdb.orig

#
# Usefull variables
#

# Get sorting correct

lipoTypes=factor(c("NF","WD","DD"),levels=c("NF","WD","DD"))

#
# Add some convience fields and fix things

crdb$TYPE=factor(crdb$MOLD_HISTOLOGY_DSCRP,levels=lipoTypes)

crdb$ASSAY="unk"
crdb$TAG=""

crdb$ASSAY[grep("^array:mRNA:Affy:U133A$",crdb$FQN)]="U133A"
crdb$ASSAY[grep("^seq:miRNA_Tuschl:Illumina:v121010$",crdb$FQN)]="miRNA"

#
# Fix filenames (used as colums) for agilent files

ii=grep("^array:aCGH:Agilent:",crdb$FQN)
crdb$ASSAY[ii]="aCGH"
crdb$TAG[ii]=gsub("(_[Ss].*de|slide|_sldie).*","",crdb$FILE_NAME[ii])
rm(ii)

####
#
# Create info for reduced set and
# remove non-Retro/Visceral samples
#

crdb.full=crdb

crdb$REDUCED.SET=NA
crdb$NUM.ASSAYS=NA

crdb=crdb[
    crdb$ASSAY!="unk"
    & !is.na(crdb$TYPE)
    & crdb$SITE1_NAME_CLEAN=="Retro-Visceral"
    & crdb$EXCLUDED!="Y",
    ]

# Fix ASSAY factor to remove 'unk' level
crdb$fASSAY=factor(crdb$ASSAY)

###################################################################################
###################################################################################
###################################################################################
#
# Get the reduced set of samples
#
# Sorting / Selection Reduced Set
#  (3) Latest MOLD_DATE (decreasing=T)
#  (2) Earliest PROCUREMENT_DATE
#  (1) Max Assay's (decreasing=T)

assayCountByBIO_ID=tapply(crdb$fASSAY,crdb$BIO_ID,function(x){len(unique(x))})
crdb$NUM.ASSAYS=assayCountByBIO_ID[match(crdb$BIO_ID,names(assayCountByBIO_ID))]

crdb=crdb[order(crdb$NUM.ASSAYS,decreasing=T),]
crdb=crdb[order(crdb$MOLD_DATE,decreasing=T),]
crdb=crdb[order(crdb$PROCUREMENT_DATE),]

uuid.assay=paste(crdb$P_MRN,crdb$TYPE,crdb$fASSAY,sep=":")
crdb=crdb[order(uuid.assay),]
uuid.assay=paste(crdb$P_MRN,crdb$TYPE,crdb$fASSAY,sep=":")
crdb$REDUCED.SET=ifelse(!duplicated(uuid.assay),"Y","")

rm(assayCountByBIO_ID)
rm(uuid.assay)

