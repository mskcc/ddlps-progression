# Methods: CGH, Cell Line, and miRNA Venn Diagram

## Figure File
`cghCellLine_miRNA_Venn_v14_.pdf`

## Generation Script
`VennTable/mkVennTable.R` (lines 349-380)

## Venn Diagram

Three-way Venn diagram showing gene assignment to each category:

- **CGH set:** Genes with positive copy number-expression concordance (CGH.RAE > 0)
- **CellLines set:** Genes with significant response to MDM2 knockdown or CDK4 inhibitor treatment, showing opposite-sign fold changes relative to tumor expression (FDR < 0.10, passed sign concordance filter)
- **miRNA set:** Genes that are predicted targets of chr12q miRNAs and show downregulation (FC < 0)

All genes are from the Set3 differentially expressed gene set. The Venn diagram displays counts for each region including individual sets and all overlaps. Detailed definitions for each category are provided in the following sections.

## Input Data

The Venn diagram displays overlaps among three gene categories derived from the integrated analysis:

### Gene Set: Set3
The analysis uses the Set3 gene set as input (defined at line 227), which consists of genes meeting the following criteria:

- Significant differential expression in U133A microarray (FDR < 0.10, |FC| > 1.3)
- Significant in at least one RNA-seq platform (pooled or paired, FDR < 0.10, |FC| > 1.3)
- Minimum of 2 platforms showing significance
- Genes ordered by sum of log2(FDR) across all FDR columns

## Three-Way Venn Diagram Categories

### Category 1: CGH
**Definition:** Genes with positive copy number-expression concordance

**Selection criteria (line 363):**
```r
CGH.RAE > 0
```

**Data source:** `CGHGenes/geneCGHvsU133aConcordence_v3.txt`

**Calculation:** From getCGHGenes.R (lines 14-23), RAE.GENE represents concordance between CGH call signs and expression fold change signs:

- If concordant samples > discordant: RAE.GENE = proportion concordant (0 to 1)
- If discordant: RAE.GENE = -1

Only genes with positive RAE.GENE values (concordant) are included in this category.

### Category 2: CellLines
**Definition:** Genes showing significant response in cell line perturbation experiments with opposite-sign concordance

**Selection criteria (line 364):**
```r
!is.na(FDR.MDM2) | !is.na(FDR.CDK4)
```

**Data sources:**

- MDM2 knockdown: `CellLines/diffGenes_20170428_MDM2_KD_380-Control_FDR_1.01_FC_1_.txt`
- CDK4 inhibitor: `CellLines/diffGenes_20170428_CDK4_inhib-Untreated_FDR_1.01_FC_1_.txt`

**Filtering applied:**

1. Initial filter: FDR < 0.10 in cell line experiment
2. Sign concordance filter: Genes must show opposite signs between tumor expression (FC.U133A) and cell line response (FC.MDM2 or FC.CDK4)
    - For MDM2 (lines 293-295): `sign(FC.U133A) != -sign(FC.MDM2)`
    - For CDK4 (lines 322-324): `sign(FC.U133A) != -sign(FC.CDK4)`

3. Genes failing sign concordance have FDR set to NA

Genes included in this category passed both significance and sign concordance filters for at least one cell line experiment.

### Category 3: miRNA
**Definition:** Genes annotated as chromosome 12q miRNA targets with consistent expression changes

**Selection criteria:**
```r
miRNA.TARGET == "X"
```

**Data source:** `Chr12q_miRNA_Targets/chr12_miRNA_ConsistentTargets.txt`

**"Consistent" definition (from `Chr12q_miRNA_Targets/get12qTargets.R`):**
Genes meeting two criteria:

1. Predicted as targets of chr12q miRNAs (hsa-mir-26a-2, hsa-mir-1279, hsa-mir-3913-1, hsa-mir-3913-2)
2. Show negative fold change (FC < 0) in U133A differential expression

The term "consistent" refers to directional consistency: genes targeted by chr12q miRNAs show downregulation as expected. The gene-miRNA target mapping is loaded from `data/db/targetsUnion.Rdata`.

**Processing:** Column "chr12.miRNA.Target.Consistent" converted to binary:

- Value of 0 → "" (blank)
- Non-zero value → "X"
