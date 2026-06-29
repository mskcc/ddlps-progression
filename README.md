# Longitudinal study of liposarcoma genomics

Crago AM, *et al*.

Code to regenerate the computational figures and tables in:

> *A multi-platform genomic analysis of paired well- and dedifferentiated
> liposarcoma identifies recurrent genomic alterations associated with
> progression and recurrence* (Crago, Lofthus, … Socci, Singer).

## Data availability

| Data | Repository | Accession |
|---|---|---|
| U133A gene expression + aCGH | GEO | **GSE244163** |
| Custom-capture validation cohort | MSKCC cBioPortal | `lipo_cbe_singers_4610_fi` |

## Processed data (Zenodo)

The precompiled snapshots and clinical export the scripts read from (`data/db/` and
`data/raw/`) are distributed separately on Zenodo: <https://doi.org/10.5281/zenodo.21041099>.

After cloning this repo, download the Zenodo tarball into the repo root and extract it
in place:

```bash
cd ddlps-progression
tar -xzf ddlps-progression-data.tar.gz
```

This will populate `data/db/` and `data/raw/` alongside the existing code.

## How to run

Scripts use relative paths and are run **from their own directory** in batch mode:

```bash
cd VennTable && Rscript --no-save mkVennTable.R
```

Many scripts rely on helper functions defined in `~/.Rprofile` (`cc()`, `write.xls()`,
`DATE()`, `len()`, etc.); source it first in a fresh session.

R 4.2.2 or newer. Required packages: `tidyverse`, `limma`, `edgeR`, `gplots`,
`IRanges`/`GenomicRanges`, `data.table`, `readxl`, `openxlsx`, `digest`, `RSQLite`,
`stringr`, `org.Hs.eg.db`, `AnnotationDbi`, `knitr`. A subset of scripts also
need `bedr` (with `bedtools` on `PATH`) and `tidygenomics`.

Analyses load processed data from `data/`: precompiled snapshots in `data/db/` and the
committed clinical export in `data/raw/CRDB/`. Each analysis directory has a
`data -> ../data` (or `../../data`) symlink so `data(...)` and `source("data/...")`
resolve.

## Figure / table → script map

| Paper item | Script | Output |
|---|---|---|
| Fig 1A heatmap | `figures/mRNAHeatmap/heatmapV2.R` | `heatMap_v14.pdf` |
| Fig 1B + Supp Table 3 | `VennTable/mkVennTable.R` | `cghCellLine_miRNA_Venn_v14_.pdf`, `joinTableCragoProgression_v14_.{txt,xlsx}` |
| Fig 2A / Supp Fig 1 (genome CNA freq) | `figures/CGHProfiles/plotRAEProfile.R` | `raePlotA0D0_{WD,DD}_*` |
| Fig 3B / 4A (6q amplification freq) | `figures/CGHProfiles/plotChrRegion.R` | `chr6Region_WD+DD_v4.pdf` |
| Fig 3A + Supp Table 7 (WD 6q+/−) | `analysis/mRNAvsCGH/integrate_cgh_mrna__WDonly__v2.R` | `suppTable_4_*_v2.xlsx` |
| Fig 6A inputs (chr13/chr8 progression blocks) | `analysis/ProgressionBlocks/findProgressionBlocks.R` | `progressionBlocks_WDvsDD_Expr_v1*.xlsx` (final image: IGV) |
| Fig 6B + Supp Tables 8/13 (DD) | `analysis/mRNAvsCGH/integrate_cgh_mrna__DDonly__v2.R` | `suppTable_12_*_v3.xlsx` |
| Supp Table 4 (per-sample CNA by region) | `tables/CGHEventTable/getCGHEventTableSelected.R` | `cghEventTable___Manuscript__SelectRegions1.xlsx` |
| Supp Table 5 (gene-level CNA in regions) | `tables/CGHEventTable/getRAEGeneTableSelect.R` | `cghGeneTable___Manuscript___SelectRegions1.xlsx` |
| Supp Tables 13/14 (band gene freq) | `tables/CGHEventTable/getRAEGeneTableByBands.R` | `raeGeneTableSigRegionsV1.xlsx` |
| JUN 1p32 freq (Supp Table 14 / text) | `analysis/JUN_Freq/cghJUN.R` | console |
| chr12q boundary (text) | `analysis/Chr12qEvent/getChr12q_Boundry.R` | `chr12qEvent.txt` |
| Supp Tables 9/10 (WES test-set mutations) | `data/maf.R` | mutation tables |
| Table 1 / Supp Tables 1,2 (clinical) | `reports/getSampleTable.R`, `reports/mkJoinTbl.R` | `sampleTableProgressionV10*.xlsx` |

## Verification

`verify.sh` at the repo root reproduces every computational output and compares it to
the committed baseline with float-tolerant diffs. See `VERIFY.md` for the step-by-step
protocol.
