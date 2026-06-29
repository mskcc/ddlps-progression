# Figure 1B — region counts and gene lists (v14, definitive)

Disjoint region counts **and the genes in each region** for the four sets in the
WDLS integrated-genomics diagram, computed from the canonical join table
`joinTableCragoProgression_v14_.txt` (produced by `VennTable/mkVennTable.R`).
Every gene is counted exactly once; all regions sum to 686.

## The four sets

| Set | Meaning | Column rule | Total |
|---|---|---|---|
| **CNA** (green) | gene in any recurrent copy-number region | `CGH.RAE > 0` | 322 |
| **12q** | the 12q13-15 amplicon — **a subset of CNA** | `CGH.Chr12 == "X"` | 20 |
| **shMDM2** | dysregulated by MDM2 knockdown (WD4847-2) | `FDR.MDM2` not NA | 72 |
| **CDK4i** | dysregulated by palbociclib (CDK4 inhibition) | `FDR.CDK4` not NA | 82 |

`12q` is wholly contained in `CNA` (all 20 12q genes have CGH.RAE>0; 0 fall outside CNA).
Cell-line union (shMDM2 OR CDK4i) = **107** = 72 + 82 - 47.

## Disjoint regions (sum = 686), with genes

Large regions (>40 genes) have their lists omitted for readability; recompute from
the join table if needed.

### Outside CNA
| Region | n | Genes |
|---|---|---|
| Other (in no set) | 301 | *(301 genes; list omitted)* |
| shMDM2 only | 16 | ABCA1, ACSL4, ACVR1, ALCAM, ANGPTL4, CYBRD1, DPP4, DSTN, FZD2, HMBS, MET, PAMR1, RGS17, RND3, SDC4, UGP2 |
| CDK4i only | 21 | ARHGAP33, C10orf10, C2CD2, CEBPD, COQ2, CTPS1, DLEU1, FBLN2, GPX3, HOOK2, IDS, IL13RA1, KLF4, MTF2, MYL9, NACC2, NUPR1, SAT1, SCPEP1, SHCBP1, YPEL5 |
| shMDM2 + CDK4i | 26 | ASPM, BUB1B, CAST, CCNB1, CCNB2, CDK1, CENPF, CKAP2, DTL, HMGB3, KIAA0101, KIF11, KIF4A, KNTC1, MCM2, MELK, NETO2, NUSAP1, ORC6, RRAGC, SMC4, TMPO, TOP2A, TYMS, XPO1, ZWINT |

### Inside CNA, but not 12q
| Region | n | Genes |
|---|---|---|
| CNA only | 262 | *(262 genes; list omitted)* |
| CNA + shMDM2 | 9 | AP3M2, BNIP3L, CDKN1A, ITGB5, LRRC32, MT1X, SDC1, SHMT2, SLC2A3 |
| CNA + CDK4i | 11 | ARID5B, ASTN2, AVPI1, CAT, IFITM2, NUP43, PCK2, PEG10, PKIA, TMBIM1, ZDHHC14 |
| CNA + shMDM2 + CDK4i | 20 | ANPEP, AURKA, BARD1, BIRC5, FANCI, FGF2, HMMR, LRP10, MAP1A, MCM4, NCAPG, PTTG1, RRM2, STMN1, TMEM97, TPX2, TRIP13, UBE2C, VAT1, VEGFB |

### Inside 12q (which is inside CNA)
| Region | n | Genes |
|---|---|---|
| 12q only | 16 | AVIL, CAND1, CCT2, CDK4, CNOT2, CTDSP2, CYP27B1, DYRK2, FRS2, KCNMB4, LGR5, MDM1, METTL1, SLC35E3, TSFM, TSPAN8 |
| 12q + shMDM2 | 0 | — |
| 12q + CDK4i | 3 | CPSF6, NUP107, YEATS4 |
| 12q + shMDM2 + CDK4i | 1 | MDM2 |

**SUM = 686**

## The 12q amplicon genes (all 20), annotated

| Gene | on shMDM2? | on CDK4i? |
|---|---|---|
| AVIL | — | — |
| CAND1 | — | — |
| CCT2 | — | — |
| CDK4 | — | — |
| CNOT2 | — | — |
| CPSF6 | — | YES |
| CTDSP2 | — | — |
| CYP27B1 | — | — |
| DYRK2 | — | — |
| FRS2 | — | — |
| KCNMB4 | — | — |
| LGR5 | — | — |
| MDM1 | — | — |
| MDM2 | YES | YES |
| METTL1 | — | — |
| NUP107 | — | YES |
| SLC35E3 | — | — |
| TSFM | — | — |
| TSPAN8 | — | — |
| YEATS4 | — | YES |

### 12q overlap with the cell-line experiments (the contested numbers)

The four 12q genes that overlap the cell-line sets, `12q * (MDM2 + CDK4)`, with
explicit per-gene membership:

| Gene | in MDM2 (shMDM2) | in CDK4 (CDK4i) |
|---|---|---|
| CPSF6 | – | YES |
| MDM2 | YES | YES |
| NUP107 | – | YES |
| YEATS4 | – | YES |

Summarised by region:

| 12q overlap | Genes | Count |
|---|---|---|
| 12q ∩ shMDM2 (total) | MDM2 | 1 |
| 12q ∩ CDK4i (total) | CPSF6, MDM2, NUP107, YEATS4 | 4 |
| 12q ∩ shMDM2 ∩ CDK4i | MDM2 | 1 |
| 12q ∩ shMDM2 only (not CDK4i) | — | 0 |
| 12q ∩ CDK4i only (not shMDM2) | CPSF6, NUP107, YEATS4 | 3 |

**Note:** only **MDM2 itself** is a 12q gene that responds to shMDM2 knockdown, and it
also responds to CDK4i (so it sits in the triple overlap). The four 12q genes that
overlap the cell-line sets — MDM2, CPSF6, YEATS4, NUP107 — are all **CDK4i** hits; only
MDM2 is also shMDM2. A figure labelling these 4 as a *shMDM2* overlap is mislabelled;
they are the **12q ∩ CDK4i** overlap. The `12q ∩ shMDM2 only` region is genuinely 0.

## Set-total cross-checks

- CNA (green, all) = **322** (of which 12q = **20**)
- shMDM2 = **72**, CDK4i = **82**, union = **107** (overlap 47)

## Notes

- These are the reproducible v14 numbers. The previously published/hand-drawn Figure 1B
  used a cell-line breakdown matching neither v14 nor v24 (see
  `RECONCILIATION_2026-06-15_v2.md` Issue 3 and `methods_fig1B.md`).
- Regenerate with `cd VennTable && Rscript --no-save mkFig1B.R`, or recompute directly
  from `joinTableCragoProgression_v14_.txt`.
- Generated 2026-06-15 from joinTableCragoProgression_v14_.txt.
