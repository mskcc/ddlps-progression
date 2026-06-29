require(data.table)

MAFFILE="db/Proj_04610_manu___SOMATIC_FACETS.vep.filtered.maf.gz"

header=readLines(MAFFILE,10)
header=header[grep("^#",header)]

mafColNames=scan("db/maf_colnames","",quiet=T)
mafColClasses=read.delim("db/maf_colClasses",header=F)
cClasses=mafColClasses$V2
names(cClasses)=mafColNames[mafColClasses$V1]

# gzip -cd (not zcat): BSD/macOS zcat insists on a .Z suffix and will not read
# a .gz file; gzip -cd reads .gz on both macOS and Linux.
maf=fread(paste("gzip -cd",MAFFILE,"|grep -v '^#'"),colClasses=cClasses)
# Why are they GL chromosomes here
maf=maf[-grep("^GL",maf$Chromosome),]

maf$BIOID=gsub("^s_","",maf$Tumor_Sample_Barcode)

maf.bioid=unique(maf$BIOID)

RM_CRDB=T
if(exists("crdb")){
    RM_CRDB=F
} else {
    source("crdb.R")
}

#
# Subset to (SITE==retro-visceral) & (! CellLines)
#

validMAFBioIDs=unique(
    crdb.full$BIO_ID[crdb.full$SITE1_NAME_CLEAN=="Retro-Visceral"
    & crdb.full$MOLD_HISTOLOGY_DSCRP %in% c("DD","WD")]
    )

reducedSamplesBIOID=maf.bioid[(maf.bioid %in% validMAFBioIDs)]
maf=maf[which(maf$BIOID %in% reducedSamplesBIOID),]

# Fix chromosome ordering by making it a factor with the correct order
maf$Chromosome=factor(maf$Chromosome,levels=c(seq(22),"X","Y","MT"))

mafs=list()

mafs$vep=list()
mafs$vep$functional=list(Consequences=c(
"frameshift_variant", "inframe_deletion", "inframe_insertion",
"missense_variant",
"regulatory_region_variant", "splice_acceptor_variant", "splice_donor_variant",
"splice_region_variant", "start_lost", "stop_gained"
))

mafs$vep$onGene=list(Variant_Classification=c(
    "3'Flank", "3'UTR", "5'Flank", "5'UTR", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Missense_Mutation", "Nonsense_Mutation", "Silent",
    "Splice_Region", "Splice_Site"))


#
# Added pseudoGermline and HQ events filter
#

maf$HQ_FILTER=""

pseudoGermline.ii=!(
    maf$Matched_Norm_Sample_Barcode=="s_FROZEN_Pooled_Normal"
    & abs(maf$t_var_freq-.5)<.05
    & maf$Existing_variation!="")

hq.ii=maf$FILTER=="."

maf$HQ_FILTER[hq.ii & pseudoGermline.ii]="X"

maf$FUNCTIONAL[maf$Consequence %in% mafs$vep$functional$Consequences]="X"

mafs$complete=maf


#
# Get patient level MAF
#

pmaf=maf
pmaf$P_MRN=crdb.full$P_MRN[match(pmaf$BIOID,crdb.full$BIO_ID)]
pmaf$PTAG=paste0(pmaf$TAG,":p_",pmaf$P_MRN)

# Sort by t_alt_count to get the "best" sample
# when we deduplicated which is the UNION method for
# collapsing over samples in a patient

pmaf=pmaf[order(t_alt_count),]
pmaf=pmaf[!duplicated(PTAG),]

mafs$patient=pmaf


if(RM_CRDB) {
    rm(list=c("CRDB.FILE", "crdb.full", "crdb.orig"))
}

rm(list=c("RM_CRDB","maf","pmaf","cClasses",
    "header", "lipoTypes", "maf.bioid", "mafColClasses", "mafColNames","crdb","pseudoGermline.ii",
    "hq.ii",
    "MAFFILE", "reducedSamplesBIOID"))
