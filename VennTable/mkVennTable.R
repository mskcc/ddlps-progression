#
# ===================================================================
# INTEGRATED GENOMIC DATA TABLE GENERATION  (VERSION 1 / v14)
# ===================================================================
#
# *** CANONICAL published version (Figure 1B + Supplemental Table 3). ***
#
# Self-contained: cell-line differential expression is read from the two
# committed CellLines/diffGenes_20170428_*.txt files. The published Fig 1B
# Venn totals come from this version's join table (cell-line union = 107,
# CDK4i = 31, 12q = 20, 686 genes total). figures/mRNAHeatmap/heatmapV2.R
# (Fig 1A) reads this script's joinTableCragoProgression_v14_.txt output.
#
# This is the only version. A later refreshed-cell-line re-run (mkVennTableV2.R,
# v24) was removed from the branch: nothing consumed its output, it was never
# the published panel, and its cell-line swap dropped MDM2 itself from the
# responder set. See RECONCILIATION_2026-06-15_v2.md.
#
# This script integrates multiple genomic data types into a unified
# gene-level table with Venn diagrams showing overlaps across platforms.
#
# Data types integrated:
# 1. Gene expression: U133A microarray, RNA-seq (pooled and paired)
# 2. Copy number: array CGH with expression concordance scores
# 3. Cell line perturbation: MDM2 knockdown, CDK4 inhibitor treatment
# 4. miRNA targets: chromosome 12q miRNA target annotations
# 5. shRNA screen: essentiality scores
#
# Output: Integrated gene table with all data types as columns, ranked
# by combined FDR scores, plus Venn diagrams showing platform overlaps.
#

VERSION="v14"

library(gplots)
library(limma)

#
# === DATA LOADING ===
# Load differential expression results, copy number data, cell line
# perturbation experiments, and functional genomics data.
#
# File naming indicates experimental comparisons:
# - U133A/RNA-seq files: "WD_vs_NF" (file names)
# - MDM2 file: "MDM2_KD_380-Control" (file name)
# - CDK4 file: "CDK4_inhib-Untreated" (file name)
# - CGH file: Contains RAE.GENE and Chr12q columns
# - miRNA file: "chr12_miRNA_ConsistentTargets" (file name)
# - shRNA: Loaded from data() function
#

u133=read.delim("U133A_NFvsWD/u133A_WDrs_vs_NF_Pfilter_ALL_v4.txt")
rnaSeqPool=read.delim("RNASeq_Progression/rnaSEQ_FIRST_WD_vs_NF_Pop__FDR_1.01_.txt")
rnaSeqPaired=read.delim("RNASeq_Progression/rnaSEQ_FIRST_WD_vs_NF_Paired__FDR__v2_1.01_.txt")
cellLine.MDM2=read.delim("CellLines/diffGenes_20170428_MDM2_KD_380-Control_FDR_1.01_FC_1_.txt")
cellLine.CDK4=read.delim("CellLines/diffGenes_20170428_CDK4_inhib-Untreated_FDR_1.01_FC_1_.txt")
cghGenes=read.delim("CGHGenes/geneCGHvsU133aConcordence_v3.txt")
data(shRNA)
miRNA=read.delim("Chr12q_miRNA_Targets/chr12_miRNA_ConsistentTargets.txt")

#
# === MASTER GENE LIST CREATION ===
# Create a union of all gene symbols measured across all platforms.
# This ensures genes detected on any platform can be evaluated for
# evidence across all other platforms.
#

genes=unique(u133$SYMBOL)
genes=union(genes,unique(rnaSeqPool$X))
genes=union(genes,unique(rnaSeqPaired$X))
genes=union(genes,unique(cghGenes$SYMBOL))
genes=union(genes,unique(cellLine.MDM2$SYMBOL))
genes=union(genes,unique(cellLine.CDK4$SYMBOL))
genes=union(genes,unique(shRNA$dd$SYMBOL))
genes=union(genes,unique(miRNA$SYMBOL))
genes=sort(genes)

#
# Initialize master integration table (vennTable) with one row per gene.
# Each platform's data will be added as columns to this table.
#

vennTable=data.frame(GENE=genes)
rownames(vennTable)=vennTable$GENE

#
# === U133A MICROARRAY DATA INTEGRATION ===
# Integrate differential expression results from U133A microarray.
#
# Multi-probe handling: Data sorted by adj.P.Val, then deduplicated by SYMBOL.
# This selects the first occurrence (lowest FDR) when multiple probes map to
# the same gene symbol.
#
# Columns added to vennTable:
# - PROBE.U133A: Probe identifier from input file column "X"
# - FDR.U133A: Adjusted p-value (adj.P.Val from input file)
# - FC.U133A: Fold change (FC from input file)
#

u133=u133[order(u133$adj.P.Val),]
u133=u133[!duplicated(u133$SYMBOL),]

vennTable$PROBE.U133A=NA
vennTable$FDR.U133A=NA
vennTable$FC.U133A=NA
ii=which(!is.na(u133$SYMBOL))
vennTable[u133$SYMBOL[ii],"PROBE.U133A"]=u133$X[ii]
vennTable[u133$SYMBOL[ii],"FDR.U133A"]=u133$adj.P.Val[ii]
vennTable[u133$SYMBOL[ii],"FC.U133A"]=u133$FC[ii]


#
# === RNA-SEQ DATA INTEGRATION ===
# Integrate RNA-seq differential expression results from two analyses:
#
# 1. POOLED samples: File name indicates "Pop" (population/pooled analysis)
#
# 2. PAIRED samples: File name indicates "Paired". The analysis script
#    diffRNASeqPAIRED.R uses design matrix ~Patient+Tissue (line 36), indicating
#    paired design with patient as blocking factor.
#
# Columns added to vennTable for each:
# - FDR.RNA_POOLED / FDR.RNA_PAIRED: False discovery rate from input files
# - FC.RNA_POOLED / FC.RNA_PAIRED: Fold change from input files
#

vennTable$FDR.RNA_POOLED=NA
vennTable$FC.RNA_POOLED=NA
vennTable[rnaSeqPool$X,"FDR.RNA_POOLED"]=rnaSeqPool$FDR
vennTable[rnaSeqPool$X,"FC.RNA_POOLED"]=rnaSeqPool$FC

vennTable$FDR.RNA_PAIRED=NA
vennTable$FC.RNA_PAIRED=NA
vennTable[rnaSeqPaired$X,"FDR.RNA_PAIRED"]=rnaSeqPaired$FDR
vennTable[rnaSeqPaired$X,"FC.RNA_PAIRED"]=rnaSeqPaired$FC

#
# === VENN DIAGRAM: GENE DETECTION ACROSS mRNA PLATFORMS ===
# Visualize which genes were measured on each expression platform.
# Converts FDR columns to binary (0 if NA, 1 if present).
#
# Gene filtering: Genes with ";" in their symbol are excluded.
#
# Output: PDF showing overlap of genes measured on U133A, RNA_POOLED, and
# RNA_PAIRED platforms.
#

#geneVenn=read.delim("geneVennTable.txt",row.names=1)
#commonGenes4=rownames(geneVenn)[rowSums(geneVenn)==4]
commonGenes=T

fdr3=vennTable[commonGenes,grep("FDR.",colnames(vennTable))]

xy=ifelse(is.na(fdr3),0,1)
xyz=xy[-grep(";",rownames(xy)),]
gsub("FDR.","",colnames(xyz))->colnames(xyz)
xyz=xyz[,c(2,3,1)]

pdf(cc("mrnaGeneSetsVenn",VERSION,".pdf"))
vennDiagram(xyz,main="mRNA Gene Sets")
dev.off()

#
# === VENN DIAGRAMS: DIFFERENTIALLY EXPRESSED GENES ===
# Identify and visualize genes showing significant differential expression
# across platforms.
#
# Significance thresholds:
# - FDR < 0.10 (qCut variable)
# - |FC| > 1.3 (FCcut variable)
#
# Creates binary matrix where 1 = significant, 0 = not significant.
# Signed matrix (sxx) multiplies by sign(FC) to track direction:
# - Positive values = upregulated and significant
# - Negative values = downregulated and significant
# - Zero = not significant
#
# Output: PDF with Venn diagrams showing gene overlaps across platforms,
# including signed version showing up/down regulation separately.
#

fc3=vennTable[commonGenes,grep("FC.",colnames(vennTable))]

fdr3[is.na(fdr3)]=1
fc3[is.na(fc3)]=1

qCut=0.10
FCcut=1.3

xx=apply((fdr3<qCut & abs(fc3)>FCcut),c(1,2),as.numeric)
sxx=xx*sign(fc3)
gsub("FC.","",colnames(sxx))->colnames(sxx)

pdf(cc("vennDiagrams_mRNA",VERSION,".pdf"))

vennDiagram(sxx[,c(2,3,1)],main="Union Genes")
vennDiagram(sxx[,c(2,3,1)],
    include=c("up","down"),
    counts.col=c("darkred","darkgreen"),
    main="Union Genes")

#
# === VENN DIAGRAMS: INTERSECTION GENES ===
# Focus on genes measured on BOTH U133A AND at least one RNA-seq platform
# (pooled or paired).
#
# Filter logic: x[1] & (x[2] | x[3]) where columns are:
# - x[1] = U133A (FDR.U133A)
# - x[2] = RNA_POOLED (FDR.RNA_POOLED)
# - x[3] = RNA_PAIRED (FDR.RNA_PAIRED)
#
# Output: Additional Venn diagrams in same PDF showing only genes measured
# across platforms.
#

fdr3=vennTable[commonGenes,grep("FDR.",colnames(vennTable))]
sixx=sxx[apply(!is.na(fdr3),1,function(x){x[1] & (x[2] | x[3])}),]
gsub("FC.","",colnames(sixx))->colnames(sixx)

vennDiagram(sixx[,c(2,3,1)],main="Intersection Genes")
vennDiagram(sixx[,c(2,3,1)],
    include=c("up","down"),
    counts.col=c("darkred","darkgreen"),
    main="Intersection Genes")

dev.off()

#
# === SET3 GENE SELECTION AND PRIORITIZATION ===
# Define "Set3" gene set from intersection genes.
#
# Selection criteria from code:
# 1. rowSums(abs(sixx)) >= 2: At least 2 platforms significant
# 2. sixx[,1] != 0: U133A must be significant (column 1 is U133A)
#
# Gene ranking: Genes ordered by sum of log2(FDR) across FDR columns.
# Lower sums (more negative log2 values) indicate stronger significance.
#

genes.set3=rownames(sixx)[rowSums(abs(sixx))>=2 & sixx[,1]!=0]

joinTbl=vennTable[genes.set3,]

# Order by sum of log(FDR)
joinTbl=joinTbl[order(rowSums(log2(joinTbl[,grep("FDR",colnames(joinTbl))]))),]

#
# === COPY NUMBER INTEGRATION: ARRAY CGH DATA ===
# Integrate copy number data from cghGenes file.
#
# RAE.GENE score: From getCGHGenes.R (lines 14-23), RAE.GENE is calculated as
# concordance between CGH call signs and expression FC signs:
# - If concordant > discordant: RAE.GENE = proportion concordant (0 to 1)
# - If discordant: RAE.GENE = -1
# Positive values indicate concordance between copy number and expression direction.
#
# Chr12q annotation: From getCGHGenes.R (lines 30-43), Chr12q column contains
# "X" if gene overlaps with chr12qEvent region, otherwise blank.
#
# Filter: Only genes with RAE.GENE > 0 are added (positive concordance).
#
# Columns added to joinTbl:
# - CGH.RAE: RAE.GENE score from input file
# - CGH.Chr12: Chr12q annotation ("X" or blank)
#

joinTbl$CGH.RAE=0
joinTbl$CGH.Chr12=""
ii.cgh=which(cghGenes$RAE.GENE>0)
mm=match(joinTbl$GENE,cghGenes$SYMBOL[ii.cgh])
nn=which(!is.na(mm))
mm=mm[nn]
joinTbl$CGH.RAE[nn]=cghGenes$RAE.GENE[ii.cgh][mm]
joinTbl$CGH.Chr12[nn]=cghGenes$Chr12q[ii.cgh][mm]

#
# === MDM2 KNOCKDOWN CELL LINE EXPERIMENT ===
# Integrate cell line perturbation data from MDM2 knockdown experiment.
#
# File name indicates: "MDM2_KD_380-Control"
#
# Sign concordance filter (lines 327-329):
# - Compares sign(FC.U133A) with sign(FC.MDM2)
# - Requires opposite signs: sign(FC.U133A) != -sign(FC.MDM2)
# - Genes with same sign have FDR.MDM2 and FC.MDM2 set to NA/0
#
# Significance filter: Only genes with FDR < qCut (0.10) are included initially.
#
# Columns added to joinTbl:
# - FDR.MDM2: False discovery rate from cell line experiment
# - FC.MDM2: Fold change from cell line experiment
# Both set to NA/0 if sign concordance filter fails.
#

joinTbl$FDR.MDM2=NA
joinTbl$FC.MDM2=0
joinTbl$FDR.CDK4=NA
joinTbl$FC.CDK4=0

ii.mdm2=which(cellLine.MDM2$FDR<qCut)
mm=match(joinTbl$GENE,cellLine.MDM2$SYMBOL[ii.mdm2])
nn=which(!is.na(mm))
mm=mm[nn]
joinTbl$FDR.MDM2[nn]=cellLine.MDM2$FDR[ii.mdm2][mm]
joinTbl$FC.MDM2[nn]=cellLine.MDM2$FC[ii.mdm2][mm]
bad.sign=which(sign(joinTbl$FC.U133A)!= -sign(joinTbl$FC.MDM2))
joinTbl$FDR.MDM2[bad.sign]=NA
joinTbl$FC.MDM2[bad.sign]=0

#
# === CDK4 INHIBITOR CELL LINE EXPERIMENT ===
# Integrate cell line perturbation data from CDK4 inhibitor treatment.
#
# File name indicates: "CDK4_inhib-Untreated"
#
# Sign concordance filter (lines 356-358):
# - Compares sign(FC.U133A) with sign(FC.CDK4)
# - Requires opposite signs: sign(FC.U133A) != -sign(FC.CDK4)
# - Genes with same sign have FDR.CDK4 and FC.CDK4 set to NA/0
#
# Significance filter: Only genes with FDR < qCut (0.10) are included initially.
#
# Columns added to joinTbl:
# - FDR.CDK4: False discovery rate from cell line experiment
# - FC.CDK4: Fold change from cell line experiment
# Both set to NA/0 if sign concordance filter fails.
#

ii.cdk4=which(cellLine.CDK4$FDR<qCut)
mm=match(joinTbl$GENE,cellLine.CDK4$SYMBOL[ii.cdk4])
nn=which(!is.na(mm))
mm=mm[nn]
joinTbl$FDR.CDK4[nn]=cellLine.CDK4$FDR[ii.cdk4][mm]
joinTbl$FC.CDK4[nn]=cellLine.CDK4$FC[ii.cdk4][mm]
bad.sign=which(sign(joinTbl$FC.U133A)!= -sign(joinTbl$FC.CDK4))
joinTbl$FDR.CDK4[bad.sign]=NA
joinTbl$FC.CDK4[bad.sign]=0

#
# === CHROMOSOME 12q miRNA TARGETS ===
# Integrate miRNA target data from chr12q miRNA targets file.
#
# File name indicates: "chr12_miRNA_ConsistentTargets"
#
# Input file contains column "chr12.miRNA.Target.Consistent" which is converted
# to binary annotation (lines 385-386):
# - Value of 0 in input -> "" (blank)
# - Non-zero value -> "X"
#
# Column added to joinTbl:
# - miRNA.TARGET: "X" if gene is chr12 miRNA target, blank otherwise
#

joinTbl$miRNA.TARGET=""

mm=match(joinTbl$GENE,miRNA$SYMBOL)
nn=which(!is.na(mm))
mm=mm[nn]
joinTbl$miRNA.TARGET[nn]=miRNA$chr12.miRNA.Target.Consistent[mm]
joinTbl$miRNA.TARGET=ifelse(joinTbl$miRNA.TARGET==0,"","X")

#
# === SECONDARY VENN DIAGRAM: CGH, CELL LINE, AND miRNA INTEGRATION ===
# Visualize overlap between three data types:
#
# 1. CGH: Genes with CGH.RAE > 0
# 2. CellLines: Genes with non-NA FDR.MDM2 OR non-NA FDR.CDK4
#    (passed both significance and sign concordance filters)
# 3. miRNA: Genes with miRNA.TARGET == "X"
#
# Output: PDF with Venn diagram showing overlap among these three categories.
#

venn2=list()

venn2$CGH=joinTbl$GENE[joinTbl$CGH.RAE>0]
venn2$CellLines=joinTbl$GENE[!is.na(joinTbl$FDR.MDM2) | !is.na(joinTbl$FDR.CDK4)]
venn2$miRNA=joinTbl$GENE[joinTbl$miRNA.TARGET=="X"]


genesU=joinTbl$GENE
vennMatrix=matrix(0,nrow=len(genesU),ncol=len(venn2))
colnames(vennMatrix)=names(venn2)
rownames(vennMatrix)=genesU

for(i in seq(len(venn2))){
    vennMatrix[venn2[[i]],i]=1
}
vennMatrix=data.frame(vennMatrix)

pdf(file=cc("cghCellLine_miRNA_Venn",VERSION,".pdf"))
vennDiagram(vennMatrix)
dev.off()

#
# === shRNA ESSENTIALITY SCREEN DATA ===
# Integrate shRNA screening data from shRNA$dd data frame.
#
# Columns added to joinTbl from shRNA$dd:
# - shRNA.RANK: RANK column from input data
# - shRNA.PCT.RANK: PCT.RANK column from input data
# - shRNA.MIN.KD: MIN.KD column from input data
# - shRNA.MAX.KD: MED.KD column from input data (note: column name mismatch
#   in code - variable named MAX.KD but populated from MED.KD)
#

joinTbl$shRNA.RANK=NA
joinTbl$shRNA.PCT.RANK=NA
joinTbl$shRNA.MIN.KD=NA
joinTbl$shRNA.MED.KD=NA
mm=match(joinTbl$GENE,shRNA$dd$SYMBOL)
nn=which(!is.na(mm))
mm=mm[nn]

joinTbl$shRNA.RANK[nn]=shRNA$dd$RANK[mm]
joinTbl$shRNA.PCT.RANK[nn]=shRNA$dd$PCT.RANK[mm]
joinTbl$shRNA.MIN.KD[nn]=shRNA$dd$MIN.KD[mm]
joinTbl$shRNA.MAX.KD[nn]=shRNA$dd$MED.KD[mm]

#
# === OUTPUT: INTEGRATED GENE TABLE ===
# Export the final integrated gene table.
#
# Table structure: Each row is a gene from Set3 (genes.set3), with columns from
# all integrated data types:
# - Gene identifier (GENE) and U133A probe ID
# - Differential expression: FDR and FC from U133A, RNA_POOLED, RNA_PAIRED
# - Copy number: CGH.RAE and CGH.Chr12
# - Cell line perturbations: FDR and FC from MDM2 and CDK4 experiments
# - miRNA targets: miRNA.TARGET
# - shRNA essentiality: RANK, PCT.RANK, MIN.KD, MAX.KD
#
# Gene order: Sorted by sum of log2(FDR) across FDR columns (line 260).
#
# Output files:
# - Excel format (.xlsx)
# - Tab-delimited format (.txt)
#

library(openxlsx)
write.xlsx(joinTbl,file=cc("joinTableCragoProgression",VERSION,".xlsx"),keepNA=FALSE,rowNames=FALSE)
write.xls(joinTbl,file=cc("joinTableCragoProgression",VERSION,".txt"),na="",row.names=F)




