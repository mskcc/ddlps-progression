#
# ===================================================================
# FIGURE 1B -- Venn membership table
# ===================================================================
#
# Reads the canonical v14 join table and writes a per-gene membership table
# marking which of the four Figure 1B sets each gene belongs to. One row per
# gene (all 686), an "X" if the gene is in that set, blank otherwise.
#
# Input : VennTable/joinTableCragoProgression_v14_.txt   (from mkVennTable.R)
# Output: VennTable/figure1B_VennDiagram_v260615.csv
#         VennTable/figure1B_VennDiagram_v260615.xlsx
#
# Set definitions (from the v14 join-table columns):
#   CNA       : CGH.RAE > 0        (gene in any recurrent copy-number region)
#   12q13-15  : CGH.Chr12 == "X"   (the 12q13-15 amplicon; a subset of CNA)
#   shMDM2    : FDR.MDM2 not NA    (dysregulated by MDM2 knockdown)
#   CDK4i     : FDR.CDK4 not NA    (dysregulated by palbociclib / CDK4 inhibition)
#
# This is a standalone read-only summary of the join table. It does not modify
# any existing file.
#

source("R/helpers.R")

suppress(library(tidyverse))

d <- read.delim("joinTableCragoProgression_v14_.txt",
                check.names=FALSE, stringsAsFactors=FALSE)

mark <- function(x) ifelse(x, "X", "")

inCNA <- d$CGH.RAE > 0
in12q <- d$CGH.Chr12 == "X"
inMDM <- !is.na(d$FDR.MDM2)
inCDK <- !is.na(d$FDR.CDK4)

# Number of sets each gene belongs to (0-4).
nSets <- inCNA + in12q + inMDM + inCDK

# Block key = the exact membership pattern, so genes in the same overlap region
# group together. Sort: (1) #sets descending, (2) block pattern, (3) gene name.
block <- paste0(as.integer(inCNA), as.integer(in12q),
                as.integer(inMDM), as.integer(inCDK))

out <- tibble(
    GENE       = d$GENE,
    CNA        = mark(inCNA),
    `12q13-15` = mark(in12q),
    shMDM2     = mark(inMDM),
    CDK4i      = mark(inCDK)
)
out <- out[order(-nSets, block, out$GENE), ]

write_csv(out, "figure1B_VennDiagram_v260615.csv")
write_xlsx(out, "figure1B_VennDiagram_v260615.xlsx")

# console sanity: per-column set sizes (should match the figure totals)
cat("wrote figure1B_VennDiagram_v260615.{csv,xlsx} :", nrow(out), "genes\n")
cat("  CNA      :", sum(out$CNA=="X"), "\n")
cat("  12q13-15 :", sum(out$`12q13-15`=="X"), "\n")
cat("  shMDM2   :", sum(out$shMDM2=="X"), "\n")
cat("  CDK4i    :", sum(out$CDK4i=="X"), "\n")
