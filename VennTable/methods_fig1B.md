# Figure 1B — auditable region counts

The published Figure 1B is a **hand-drawn** 4-region diagram that lives outside
this repo. Its region totals were assigned by hand and **do not close**: the
visible numbers (260, 20, 1, 11, 16, 12, 22, 19, 23, 304) sum to **688**, not
the **686** genes in the join table. `mkFig1B.R` recomputes every region from
the canonical v14 join table so the counts are reproducible and sum to 686
exactly (enforced by a `stopifnot`). It does not redraw the figure.

## Script

`VennTable/mkFig1B.R`
- reads `joinTableCragoProgression_v14_.txt` (the canonical v14 join table from
  `mkVennTable.R`),
- defines the four sets (below),
- writes `fig1B_counts_v14_.txt` (every disjoint region, auditable, sums to 686).

Run: `cd VennTable && Rscript --no-save mkFig1B.R`

## The four sets (all within "DE genes in WDLS vs normal fat")

| Set (figure color) | Meaning | Column rule in join table | Total |
|---|---|---|---|
| 12q13-15 (dark green) | gene encoded on the 12q13-15 amplicon | `CGH.Chr12 == "X"` | 20 |
| CNA-other (green) | CNA-concordant, region **other** than 12q | `CGH.RAE > 0 & not 12q` | 302 |
| shMDM2 (blue) | altered by MDM2 knockdown (WD4847-2) | `FDR.MDM2` not NA | 72 |
| CDK4i (purple) | altered by palbociclib (CDK4 inhibition) | `FDR.CDK4` not NA | 82 |
| Other (outer) | none of the four | — | 301 |

cell-line union (shMDM2 OR CDK4i) = 107.

## Published figure vs reproducible counts

The figure's intent matches these sets, but its hand numbers differ slightly
(and sum to 688). The reproducible counts (sum to 686):

| Figure region (intent) | Figure (hand) | Reproducible (v14) |
|---|---|---|
| Other (no set) | 304 | **301** |
| 12q13-15 total | 21 | **20** |
| CNA-other (green only) | 260 | **262** |
| cell-line union | 107 | **107** |
| grand total | 688 (does not close) | **686** |

Use the reproducible numbers; the prose has already been updated to 686 / 20 / 107 /
253 / 17.9%.

## Bookkeeping notes

- The `12q13-15 & shMDM2` region is genuinely 0 in the v14 data (only MDM2 itself
  is a 12q gene responding to shMDM2 knockdown, and it also responds to CDK4i, so
  it sits in the triple overlap). The hand-drawn figure conflates this with the
  `12q & CDK4i` overlap.
- The published panel nested shMDM2 *inside* a combined "shMDM2 or CDK4i" set and
  used hand counts. This script instead treats shMDM2 and CDK4i as the two real
  table columns -- the honest, reproducible decomposition. The biological message
  is identical; the bookkeeping now closes.
