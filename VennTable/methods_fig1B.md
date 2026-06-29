# Figure 1B — making it reproducible

The published Figure 1B is a **hand-drawn** nested 4-region diagram. Its region
totals were assigned by hand and **do not close**: the visible numbers
(260, 20, 1, 11, 16, 12, 22, 19, 23, 304) sum to **688**, not the **686** genes in
the join table. `mkFig1B.R` regenerates the same diagram from the canonical v14
join table so every region is reproducible and the counts sum to 686 exactly
(enforced by a `stopifnot`).

## Script

`VennTable/mkFig1B.R`
- reads `joinTableCragoProgression_v14_.txt` (the canonical v14 join table from
  `mkVennTable.R`),
- defines the four sets (below),
- writes `fig1B_counts_v14_.txt` (every disjoint region, auditable, sums to 686),
- draws an area-proportional Euler diagram with `eulerr` → `fig1B_v14_.pdf`.

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
253 / 17.9% (see `RECONCILIATION_2026-06-15_v2.md` Issue 4).

## Fidelity caveats (read before "fixing" the diagram)

- **Counts are exact; areas are best-fit.** `eulerr` solves a 4-set area-proportional
  layout numerically and cannot always render every region with perfect area. One
  pinch-point region (`12q13-15 & shMDM2` only) is genuinely 0 in the data and shows
  as "0" at a boundary. This is correct, not a bug.
- The published panel nested shMDM2 *inside* a combined "shMDM2 or CDK4i" set and
  used hand counts. `mkFig1B.R` instead treats shMDM2 and CDK4i as the two real
  table columns (the honest, reproducible decomposition). The biological message is
  identical; the bookkeeping now closes.
- This is a *figure regenerated from data*, not a pixel match to the hand drawing.
  If the journal requires the exact published styling, use this script's counts to
  annotate the existing artwork rather than replacing it.
