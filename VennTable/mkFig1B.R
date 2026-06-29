#
# ===================================================================
# FIGURE 1B  (reproducible Euler/Venn of the integrated WDLS gene set)
# ===================================================================
#
# The published Figure 1B is a hand-drawn nested 4-region diagram. Its region
# totals were assembled by hand and do not quite close (they sum to 688, not the
# 686 genes in the join table). This script regenerates the SAME diagram from the
# canonical v14 join table so the numbers are reproducible and sum to 686 exactly.
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
# Outputs (in VennTable/):
#   fig1B_v14_.pdf            -- the diagram
#   fig1B_counts_v14_.txt     -- every region count (auditable; sums to 686)
#

suppress <- suppressPackageStartupMessages
suppress(library(tidyverse))
suppress(library(eulerr))

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
# eulerr/euler counts each unique combination of memberships exactly once.
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

#
# === DIAGRAM ===
# Area-proportional Euler diagram of the four sets, drawn inside a labeled box
# whose caption carries the outer "Other" count (the genes in none of the sets),
# matching the published Figure 1B layout.
#
fit <- euler(mm)

pdf(cc("fig1B",VERSION,".pdf"), width=8, height=7)
print(plot(fit,
    quantities = list(type="counts", cex=0.9),
    fills = list(fill=c("#1b9e77","#66c2a5","#8da0cb","#9e6ebd"), alpha=0.55),
    edges = list(col="grey30"),
    labels = list(cex=0.95),
    main = sprintf(
      "Figure 1B: %d DE genes in WDLS vs normal fat\n(%d not in any set = 'Other genes with altered expression')",
      N, nOther)))
dev.off()

cat("Wrote", cc("fig1B",VERSION,".pdf"), "and", cc("fig1B_counts",VERSION,".txt"), "\n")
