#
# ===================================================================
# FIGURE 1B  (auditable region counts for the integrated WDLS gene set)
# ===================================================================
#
# The published Figure 1B is a hand-drawn 4-region diagram that lives outside
# this repo. Its region totals were assembled by hand and do not close: they
# sum to 688, not the 686 genes in the join table. This script does NOT redraw
# the figure -- it produces the auditable region-count table that exposed the
# 688 vs 686 discrepancy and that the prose now cites.
#
# Source of truth: VennTable/joinTableCragoProgression_v14_.txt (from mkVennTable.R).
# See VennTable/methods_fig1B.md for the full mapping of figure regions to columns.
#
# The four sets (all within "genes differentially expressed in WDLS vs normal fat"):
#   12q13-15   : gene encoded on the 12q13-15 amplicon      (CGH.Chr12 == "X")
#   CNA-other  : CNA-concordant, on a region OTHER than 12q  (CGH.RAE>0 & not 12q)
#   shMDM2     : altered by MDM2 knockdown in WD4847-2        (FDR.MDM2 not NA)
#   CDK4i      : altered by palbociclib (CDK4 inhibition)     (FDR.CDK4 not NA)
# Genes in none of the four = "Other genes with altered expression" (outer region).
#
# Output (in VennTable/):
#   fig1B_counts_v14_.txt     -- every region count (auditable; sums to 686)
#

suppress <- suppressPackageStartupMessages
suppress(library(tidyverse))

VERSION <- "v14"

d <- read.delim("joinTableCragoProgression_v14_.txt",
                check.names=FALSE, stringsAsFactors=FALSE)
N <- nrow(d)

#
# === SET DEFINITIONS ===
#
on12q    <- d$CGH.Chr12 == "X"
cnaOther <- d$CGH.RAE > 0 & !on12q
mdm2     <- !is.na(d$FDR.MDM2)
cdk4     <- !is.na(d$FDR.CDK4)

sets <- c("12q13-15","CNA-other","shMDM2","CDK4i")
mm <- cbind(`12q13-15`=on12q, `CNA-other`=cnaOther, shMDM2=mdm2, CDK4i=cdk4)
mode(mm) <- "logical"

inAny <- rowSums(mm) > 0
nOther <- sum(!inAny)

#
# === AUDITABLE REGION COUNTS (must sum to N=686) ===
# Each gene contributes to exactly one disjoint region (its combination of
# memberships), so region counts sum to N.
#
combo <- apply(mm, 1, function(r) {
    s <- sets[r]
    if (len(s)==0) "Other" else paste(s, collapse="&")
})
counts <- sort(table(combo), decreasing=TRUE)

con <- file(cc("fig1B_counts",VERSION,".txt"),"w")
writeLines(c(
  "Figure 1B -- reproducible region counts from joinTableCragoProgression_v14_.txt",
  sprintf("Total genes (N) = %d", N),
  "",
  "Set definitions:",
  "  12q13-15  : CGH.Chr12 == 'X'",
  "  CNA-other : CGH.RAE > 0 AND not 12q13-15",
  "  shMDM2    : FDR.MDM2 not NA",
  "  CDK4i     : FDR.CDK4 not NA",
  "  Other     : none of the above",
  "",
  "Set totals:",
  sprintf("  12q13-15        : %d", sum(on12q)),
  sprintf("  CNA-other       : %d", sum(cnaOther)),
  sprintf("  shMDM2          : %d", sum(mdm2)),
  sprintf("  CDK4i           : %d", sum(cdk4)),
  sprintf("  cell-line union : %d", sum(mdm2|cdk4)),
  sprintf("  Other           : %d", nOther),
  "",
  "Disjoint region counts (each gene counted once; these sum to N):"),con)
for (nm in names(counts))
    writeLines(sprintf("  %-30s %d", nm, counts[[nm]]), con)
writeLines(c("", sprintf("SUM = %d  (must equal %d)", sum(counts), N)), con)
close(con)

cat(readLines(cc("fig1B_counts",VERSION,".txt")), sep="\n"); cat("\n")
stopifnot(sum(counts) == N)

cat("Wrote", cc("fig1B_counts",VERSION,".txt"), "\n")
