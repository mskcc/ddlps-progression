# VERIFY.md — verifying the public `master` branch on the HPC server

This document is a step-by-step protocol to confirm that the pruned, paper-only public repo
(`master`) is intact and reproduces the manuscript's computational outputs **on the
HPC server**, where the full environment exists (R with `openxlsx`, `bedr`/bedtools,
the `~/.Rprofile` helpers, and the institutional `/ifs/rtsia01/...` paths). All Excel
writes in the figure pipeline have been migrated off the rJava-backed `xlsx` package to
`openxlsx`, so neither `xlsx` nor a Java install is required.

(Historical: the pruning work happened on branch `manu/v_2026`, which is now the `master`
branch of the public repo.)

It was written off-server, so several scripts could only be *partially* verified locally
(no `bedr`, no `/ifs` paths). Steps below are flagged:

- **[verified-local]** — already run off-server; expected to "just work" here.
- **[hpc-first-run]** — first executed here; give these extra scrutiny.

## Conventions used in every check

- Every script is run **from its own directory** with `Rscript --no-save`. The scripts use
  relative paths; running from the wrong directory will fail.
- Data is loaded through R's `data()` mechanism (explained in Step 3). The `data/`
  symlinks must be intact.
- **Reproduction is checked with git**, not by eyeballing. The committed copy of each
  output *is* the paper baseline. After running a script you compare the regenerated file
  against the committed one:
  - text output → `git diff --no-color <file>` shows content changes; `verify.sh` does a
    float-tolerant comparison (cell-by-cell, rel-tol `1e-9`), so trailing-digit float
    rounding (HPC vs Mac IEEE-754) is counted as PASS, real value changes as WARN/FAIL;
  - binary output (`.xlsx`/`.pdf`/`.rda`) → `git diff --stat <file>` shows whether bytes
    changed at all.
  - Then discard the regenerated file with `git checkout -- <file>` so the branch stays
    clean.
- **PASS / FAIL** criteria are stated per step. "Trivial diff" = only floating-point last
  digits (e.g. `79.36168315949` vs `79.3616831594901`) or PDF/xlsx timestamps. Those are
  expected and are a PASS. A real diff = changed gene lists, row counts, or values.

Set this once so the snippets are copy-pasteable:

```bash
export REPO=/path/to/ddlps-progression      # <-- edit to the clone path on HPC
cd "$REPO"
```

> **Shortcut:** `verify.sh` at the repo root automates this entire protocol.
> After loading your modules (Step 0), run `./verify.sh` to execute every step
> and print a PASS/FAIL/WARN summary table; `./verify.sh --list` shows the step
> ids; `./verify.sh 4a 6a` runs only selected steps. It restores every
> regenerated output afterward so the working tree stays clean. The manual
> steps below remain the authoritative reference for the PASS/FAIL criteria and
> for inspecting any step the script flags **WARN**.

---

## Step 0 — Environment / modules

The repo does not load HPC modules for you. Load your site's modules first, e.g.:

```bash
module load R/4.2.2        # floor is R 4.2.2; development was on 4.5.1
module load bedtools       # required for the bedr package
```

(No Java/`xlsx` module is needed — all `.xlsx` writes now go through `openxlsx`.)

PASS: the next command prints R 4.2.2 or newer. (`verify.sh` enforces the floor
numerically; override with `R_MIN=4.3.0 ./verify.sh` if needed.)

```bash
Rscript --version
```

---

## Step 1 — Branch, working tree, and `~/.Rprofile`

```bash
cd "$REPO"
git status -sb | head -1            # expect: ## master
git log --oneline -1               # the prune/manifest commit (or note if not yet committed)
ls -l ~/.Rprofile                  # must exist; defines cc(), write.xls(), DATE(), len(), etc.
```

PASS:
- branch is `master`, **or** any feature branch descended from it (`verify.sh`
  accepts both — useful when working on a `fix/*` branch off `master`);
- `~/.Rprofile` exists. The HPC `~/.Rprofile` must define `cc`, `write.xls`, `write_xlsx`,
  `DATE`, `len`, `getSDIR`, `halt`, `suppress`. Confirm:

```bash
Rscript -e 'source("~/.Rprofile"); cat(all(sapply(c("cc","write.xls","write_xlsx","DATE","len","getSDIR","halt","suppress"), exists)), "\n")'
```

PASS: prints `TRUE`. FAIL: if `FALSE`, the HPC `.Rprofile` differs from the one the code
was written against — stop and reconcile before continuing.

---

## Step 2 — Required R packages

```bash
Rscript -e '
pkgs <- c("tidyverse","limma","edgeR","gplots","IRanges","GenomicRanges",
          "data.table","readxl","openxlsx","digest","RSQLite","stringr",
          "org.Hs.eg.db","AnnotationDbi","knitr",
          "bedr","tidygenomics")
inst <- rownames(installed.packages())
for (p in pkgs) cat(sprintf("%-14s %s\n", p, if (p %in% inst) "OK" else "*** MISSING ***"))
'
```

PASS: all OK. Note specifically:
- **`openxlsx`** is the workbook writer used throughout the figure pipeline (no Java
  dependency). The migration from `xlsx` to `openxlsx` covers
  `VennTable/mkVennTable.R`, `analysis/ProgressionBlocks/findProgressionBlocks.R`,
  `tables/CGHEventTable/getCGHEventTableSelected.R`, `reports/getSampleTable.R`, and
  `reports/mkJoinTbl.R`. The only remaining `xlsx` caller in the repo is
  `data/raw/CRDB/updateCRDBWithSiteInfo.R`, a data-prep script outside the figure
  pipeline and outside `verify.sh`.
- **`bedr`** (needs `bedtools` on PATH) is required by
  `tables/CGHEventTable/getRAEGeneTableByBands.R`. Confirm the binary too:
- **`tidygenomics`** is required by `tables/CGHEventTable/getRAEGeneTableSelect.R`
  (`genome_intersect()`).

```bash
which bedtools && bedtools --version
```

---

## Step 3 — Data-loading mechanism (the foundation)

**How loading works (important):** the scripts call `data(u133a)`, `data(sampleTable)`,
etc. R's `data()` searches a `data/` subdirectory of the *current working directory*. In
this repo every analysis folder has a `data -> ../data` (or `../../data`) symlink, and the
real `data/` directory contains BOTH:
- `.rda`/`.Rdata` objects that `data()` loads directly: **`u133a.rda`** (~27 MB),
  `cellLines.rda`, `targetsUnion.Rdata`; and
- `.R` loader scripts that `data()` **executes**: `sampleTable.R`, `crdb.R`, `cghUBM.R`,
  `cghGeneMatrix.R`, `cghEventsPerBlock.R`, `rnaSeq.R`, `miRNATargets.R`,
  `progressionSet.R`, `regions.R`, `hg18.R`, `shRNA.R`, plus helpers `tools.R`, `raeCNV.R`,
  `maf.R`. These in turn read from `data/db/` (the ~180 MB precompiled snapshots) and the
  committed clinical export `data/raw/CRDB/FullDataExport_...2017_05_03.txt`.

### 3a. All `data` symlinks resolve (the prune must not have broken any)

```bash
cd "$REPO"
broken=0
for l in $(find . -type l -name data -not -path './.git/*'); do
  if [ ! -e "$l/" ]; then echo "BROKEN: $l -> $(readlink "$l")"; broken=1; fi
done
[ $broken -eq 0 ] && echo "ALL DATA SYMLINKS OK"
```

PASS: prints `ALL DATA SYMLINKS OK`.

### 3b. Required input blobs are present

```bash
cd "$REPO"
for f in \
  data/u133a.rda data/cellLines.rda data/targetsUnion.Rdata \
  data/sampleTable.csv data/chr12qEvent.csv \
  data/raw/CRDB/FullDataExport_Nick_Socci_04_05_2016___PaperFREEZE_2017_05_03.txt \
  data/raw/CGH/cghCBSRdataFiles \
  data/raw/CGH/rae/FEAT.file \
  data/db/XX_lumiNorm_20160516_.Rdata \
  data/db/CragoProgressionC-lesions.Rdata \
  data/db/GeneMatrix.txt \
  data/db/CragoProgressionC_EventMatrixV2_20160419_.txt \
  data/db/Proj_04610_manu___SOMATIC_FACETS.vep.filtered.maf.gz \
  data/db/Proj_3704_Merge_GeneCounts.txt.gz \
  data/db/targetsUnion.Rdata data/db/chr12q_mirna.txt \
  data/db/maf_colClasses data/db/maf_colnames \
  data/db/progressionEventTable.csv \
  data/db/averageSignal_3Sort_M3_.txt data/db/human.hg18.genome data/db/survialRegions \
  data/db/progressionSet.txt data/db/CNV.file \
  analysis/mRNAvsCGH/Rlib/annotation.R analysis/mRNAvsCGH/Rlib/hgu133a.sqlite ; do
  [ -e "$f" ] && echo "OK   $f" || echo "MISSING $f"
done
```

PASS: every line is `OK`. Two of these — `analysis/mRNAvsCGH/Rlib/{annotation.R,hgu133a.sqlite}`
— were **moved here from `Pass1/`** during pruning (see Step 6); their presence confirms
the move survived.

### 3c. Smoke-test the loaders **[verified-local]**

```bash
cd "$REPO/data"
Rscript --no-save -e '
source("~/.Rprofile")
suppressWarnings(suppressMessages({
  source("sampleTable.R"); source("progressionSet.R"); source("crdb.R")
  source("cghUBM.R");      source("cghGeneMatrix.R");   source("rnaSeq.R")
  source("miRNATargets.R");source("shRNA.R");           source("hg18.R")
  source("regions.R")
}))
stopifnot(nrow(sampleTable) == 349)
cat("sampleTable rows:", nrow(sampleTable), "\n")
cat("ubm$a0 dims     :", paste(dim(ubm$a0), collapse=" x "), "\n")
cat("geneCGH dims    :", paste(dim(geneCGH), collapse=" x "), "\n")
cat("rnaSeq$ds dims  :", paste(dim(rnaSeq$ds), collapse=" x "), "\n")
cat("LOADERS_OK\n")
'
```

PASS (expected numbers, verified off-server):
- `sampleTable rows: 349`
- `ubm$a0 dims : 19933 x 187`
- `geneCGH dims : 24904 x 187`
- `rnaSeq$ds dims : 20032 x 24`
- final line `LOADERS_OK`

FAIL: a different row/col count means the `data/db/` snapshots on HPC differ from those in
the branch — investigate which file.

### 3d. `data/maf.R` loads the WES MAF **[verified-local]** (Supp Tables 9/10)

> Reads the gzipped MAF via `gzip -cd` (was `zcat`, which fails on macOS BSD).
> Verified locally: `mafs$complete` 71319 rows, `mafs$patient` 60895 rows.

```bash
cd "$REPO/data"
Rscript --no-save -e '
source("~/.Rprofile")
suppressWarnings(suppressMessages(source("maf.R")))
cat("mafs$complete rows:", nrow(mafs$complete), "\n")
cat("mafs$patient  rows:", nrow(mafs$patient), "\n")
cat("MAF_OK\n")
'
```

PASS: prints two positive row counts and `MAF_OK` with no error. (Record the counts; they
become the reference for any later mutation-table work.)

---

## Step 4 — Figure 1 chain: Venn (1B) → heatmap (1A)

### 4a. `VennTable/mkVennTable.R` (v14) — Fig 1B + Supp Table 3 **[verified-local]**

The intermediate inputs it reads (`U133A_NFvsWD/...v4.txt`, `RNASeq_Progression/...txt`,
`CGHGenes/...v3.txt`, `Chr12q_miRNA_Targets/...txt`) are committed, so this runs directly.
The script writes the three Venn PDFs, the join-table `.txt`, and the workbook `.xlsx`
all through `openxlsx` (no Java required).

```bash
cd "$REPO/VennTable"
Rscript --no-save mkVennTable.R 2>&1 | tail -5
echo "---- compare regenerated outputs to committed baseline ----"
git --no-pager diff --stat -- joinTableCragoProgression_v14_.txt \
    cghCellLine_miRNA_Venn_v14_.pdf mrnaGeneSetsVenn_v14_.pdf \
    vennDiagrams_mRNA_v14_.pdf joinTableCragoProgression_v14_.xlsx
git --no-pager diff --no-color -- joinTableCragoProgression_v14_.txt | head -40
```

PASS:
- script runs to completion and writes the `.txt`, `.xlsx`, and three Venn PDFs;
- `git diff --no-color joinTableCragoProgression_v14_.txt` shows **either nothing or only
  trailing-digit float differences**. `verify.sh` does this comparison cell-by-cell with
  rel-tol `1e-9`, so HPC-vs-Mac IEEE-754 rounding noise is automatically classified PASS.

FAIL: any gene added/removed or any value differing beyond the ~13th significant digit.

Reset the regenerated files:

```bash
git checkout -- joinTableCragoProgression_v14_.txt cghCellLine_miRNA_Venn_v14_.pdf \
  mrnaGeneSetsVenn_v14_.pdf vennDiagrams_mRNA_v14_.pdf joinTableCragoProgression_v14_.xlsx 2>/dev/null
```

### 4b. Upstream DE regenerators (not run by `verify.sh`)

The seven inputs `mkVennTable.R` reads are committed to the repo, and step 4a treats them
as the boundary of the figure-reproduction check. The scripts that originally produced
those inputs are preserved alongside each input as the canonical record of its derivation,
but `verify.sh` does NOT run them:

- `VennTable/U133A_NFvsWD/doU133A_WDvsNF.R` → `u133A_WDrs_vs_NF_Pfilter_ALL_v4.txt`
- `VennTable/RNASeq_Progression/diffRNASeqPOOLED.R` → `rnaSEQ_FIRST_WD_vs_NF_Pop__FDR_1.01_.txt`
- `VennTable/RNASeq_Progression/diffRNASeqPAIRED.R` → `rnaSEQ_FIRST_WD_vs_NF_Paired__FDR__v2_1.01_.txt`
- `VennTable/CGHGenes/getCGHGenes.R` → `geneCGHvsU133aConcordence_v3.txt`
- `VennTable/Chr12q_miRNA_Targets/get12qTargets.R` → `chr12_miRNA_ConsistentTargets.txt`
  (reads the frozen sibling `../U133A_NFvsWD/u133A_WDrs_vs_NF_Pfilter_ALL_v3.txt`; no
  script regenerates `_v3.txt`)
- `VennTable/CellLines/diffCellLines.R` → `CellLines/diffGenes_20170428_*.txt`
  (the script's output filename embeds `DATE()`, so re-running it today would not produce
  the dated files `mkVennTable.R` reads — see `mkVennTable.R` lines 8-10)

To exercise any of these regenerators by hand, run it from its own directory and compare
the regenerated output against the committed copy with `git diff` (then
`git checkout --` to discard the regenerated file).

### 4c. `figures/mRNAHeatmap/heatmapV2.R` — Fig 1A **[verified-local]**

Reads `joinTableCragoProgression_v14_.txt` via a symlink to `VennTable/`.

```bash
cd "$REPO/figures/mRNAHeatmap"
ls -l joinTableCragoProgression_v14_.txt        # must be a symlink that resolves
Rscript --no-save heatmapV2.R 2>&1 | tail -3
git --no-pager diff --stat -- heatMap_v14.pdf
git checkout -- heatMap_v14.pdf
```

PASS: exits cleanly and writes `heatMap_v14.pdf` (~0.5–0.6 MB). Binary diff may show a byte
change (PDF timestamp) — open it and confirm it is the WDLS-vs-normal-fat clustered heatmap.

---

## Step 5 — Copy-number figures and tables

### 5a. `figures/CGHProfiles/plotRAEProfile.R` — Fig 2A / Supp Fig 1 **[hpc-first-run]**

Genome-wide amplification/deletion frequency profiles. Takes a sample-type argument in the
original workflow; run it as the committed outputs imply (WD and DD):

```bash
cd "$REPO/figures/CGHProfiles"
head -20 plotRAEProfile.R           # confirm how WD/DD is selected (arg vs edit)
Rscript --no-save plotRAEProfile.R 2>&1 | tail -5
ls -l raePlotA0D0_*.{pdf,png} 2>/dev/null
```

PASS: produces `raePlotA0D0_WD_*` and `raePlotA0D0_DD_*` without error. Compare against the
committed versions (`git diff --stat`), then `git checkout --` them.

### 5b. `figures/CGHProfiles/plotChrRegion.R` — Fig 3B / 4A **[verified-local]**

```bash
cd "$REPO/figures/CGHProfiles"
Rscript --no-save plotChrRegion.R 2>&1 | tail -3
git --no-pager diff --stat -- "chr6Region_WD+DD_v4.pdf" chr6RegionTest_v3.pdf chr6RegionTest_WDOnly_v3.pdf
git checkout -- "chr6Region_WD+DD_v4.pdf" chr6RegionTest_v3.pdf chr6RegionTest_WDOnly_v3.pdf
```

PASS: exits cleanly, writes `chr6Region_WD+DD_v4.pdf` (~13 KB). It is the 6q amplification-
frequency profile (WD vs DD) with SASH1/CCDC28A/TAB2 loci.

### 5c. `tables/CGHEventTable/getCGHEventTableSelected.R` — Supp Table 4 **[verified-local]**

```bash
cd "$REPO/tables/CGHEventTable"
Rscript --no-save getCGHEventTableSelected.R 2>&1 | tail -3
git --no-pager diff --stat -- cghEventTable___Manuscript__SelectRegions1.xlsx
git checkout -- cghEventTable___Manuscript__SelectRegions1.xlsx
```

PASS: writes `cghEventTable___Manuscript__SelectRegions1.xlsx` with no error.

### 5d. `tables/CGHEventTable/getRAEGeneTableSelect.R` — Supp Table 5 **[hpc-first-run]**

```bash
cd "$REPO/tables/CGHEventTable"
Rscript --no-save getRAEGeneTableSelect.R 2>&1 | tail -3
git --no-pager diff --stat -- cghGeneTable___Manuscript___SelectRegions1.xlsx
git checkout -- cghGeneTable___Manuscript___SelectRegions1.xlsx
```

PASS: writes `cghGeneTable___Manuscript___SelectRegions1.xlsx` with no error. (Reads the
committed `ucsc_hg18__*.txt.gz` reference files in this directory.)

### 5e. `tables/CGHEventTable/getRAEGeneTableByBands.R` — Supp Tables 13/14 **[hpc-first-run — needs `bedr`]**

This one failed off-server **only** because `bedr`/`bedtools` was unavailable
(`could not find function "bedr.sort.region"`). Its data dependency — the preserved
`data/raw/CGH/rae/FEAT.file` — read correctly. With `bedtools` loaded (Step 0/2) it should
complete.

```bash
cd "$REPO/tables/CGHEventTable"
Rscript --no-save getRAEGeneTableByBands.R 2>&1 | tail -6
git --no-pager diff --stat -- raeGeneTableSigRegionsV1.xlsx
git checkout -- raeGeneTableSigRegionsV1.xlsx
```

PASS: completes and writes `raeGeneTableSigRegionsV1.xlsx`. FAIL: if it still errors on
`bedr.sort.region`, the `bedtools` binary is not on PATH — fix the module, do not change
the code.

---

## Step 6 — Integration scripts (Fig 3A / 6B) — the moved `Rlib` dependency

During pruning, `Pass1/OldCopyNumPaper/Rlib/{annotation.R,hgu133a.sqlite}` were **moved**
to `analysis/mRNAvsCGH/Rlib/`, and these two scripts were repointed from
`source("../../Pass1/OldCopyNumPaper/Rlib/annotation.R")` to `source("Rlib/annotation.R")`
(and likewise for the sqlite path). These steps confirm the move + repoint are correct.
Both reproduced **identically** off-server. They use `openxlsx` (not `xlsx`).

> Note: `Rlib/annotation.R` also defines a `load.GO.ann()` that references
> `~/Work/AnnoteDBs/GO.sqlite`. The integration scripts call only `getAnnoteTable()`, **not**
> `load.GO.ann()`, so that GO path is not needed. If a future edit calls `load.GO.ann()`,
> that file would have to exist on HPC.

### 6a. WD integration — Fig 3A + Supp Table 7 **[verified-local, identical]**

```bash
cd "$REPO/analysis/mRNAvsCGH"
ls -l Rlib/annotation.R Rlib/hgu133a.sqlite survialRegions   # deps must resolve
grep -n 'Rlib/' integrate_cgh_mrna__WDonly__v2.R             # paths must say Rlib/ (not Pass1/)
Rscript --no-save integrate_cgh_mrna__WDonly__v2.R 2>&1 | tail -5
git --no-pager diff --stat -- suppTable_4_with_U133A__WDpEvent_vs_WD_without_v2.xlsx \
    diff_mRNA_vs_CopyNumber_by_blocks___WDonly_v7_.rda
git checkout -- suppTable_4_with_U133A__WDpEvent_vs_WD_without_v2.xlsx \
    diff_mRNA_vs_CopyNumber_by_blocks___WDonly_v7_.rda 2>/dev/null
```

PASS:
- `grep` shows `source("Rlib/annotation.R")` and `getAnnoteTable("Rlib/hgu133a.sqlite")`
  (NOT any `Pass1/...` path);
- script completes (benign many-to-many `left_join` warnings are expected, not errors);
- writes `suppTable_4_..._v2.xlsx`. Optional stronger check (xlsx content equality), run
  *before* the `git checkout`:

```bash
Rscript -e 'a=readxl::read_xlsx("suppTable_4_with_U133A__WDpEvent_vs_WD_without_v2.xlsx"); cat("dim:",dim(a),"\n")'
# off-server reference: dim 127 7
```

FAIL: a `cannot open file 'Rlib/...'` or SQLite error means the move/repoint is wrong.

### 6b. DD integration — Fig 6B + Supp Tables 8/13 **[verified-local]**

```bash
cd "$REPO/analysis/mRNAvsCGH"
grep -n 'Rlib/' integrate_cgh_mrna__DDonly__v2.R
Rscript --no-save integrate_cgh_mrna__DDonly__v2.R 2>&1 | tail -5
git --no-pager diff --stat -- suppTable_12_with_U133A__DDpEvent_vs_DD_without_v3.xlsx \
    diff_mRNA_vs_CopyNumber_by_blocks___DDonly_v8_.rda
git checkout -- suppTable_12_with_U133A__DDpEvent_vs_DD_without_v3.xlsx \
    diff_mRNA_vs_CopyNumber_by_blocks___DDonly_v8_.rda 2>/dev/null
```

PASS: completes; writes `suppTable_12_..._v3.xlsx`; paths reference `Rlib/`.

---

## Step 7 — Remaining analyses

### 7a. `analysis/ProgressionBlocks/findProgressionBlocks.R` — Fig 6A inputs **[verified-local]**

```bash
cd "$REPO/analysis/ProgressionBlocks"
Rscript --no-save findProgressionBlocks.R 2>&1 | tail -5
ls -l progressionBlocks_WDvsDD_Expr_v1*.xlsx 2>/dev/null
git --no-pager diff --stat -- progressionBlocks_WDvsDD_Expr_v1*.xlsx 2>/dev/null
git checkout -- progressionBlocks_WDvsDD_Expr_v1*.xlsx 2>/dev/null
```

PASS: completes and writes the progression-blocks workbook (chr13 / chr8 sheets). The final
IGV image (Fig 6A) is produced manually from these calls, not by this script.

### 7b. `analysis/JUN_Freq/cghJUN.R` — JUN 1p32 frequency **[verified-local]**

> This file was previously a raw REPL transcript (a saved interactive session with a
> half-typed broken line) that could not run with `Rscript`. It has been rewritten as a
> clean runnable script containing only the actual computation; the frequencies it prints
> are unchanged.

```bash
cd "$REPO/analysis/JUN_Freq"
Rscript --no-save cghJUN.R 2>&1 | tail -15
```

PASS: runs and prints JUN-amplification frequencies by TYPE consistent with the paper
(WDLS amplified 0.057, DDLS amplified 0.269; Supp Table 14). Console/table script —
capture the output.

### 7c. `analysis/Chr12qEvent/getChr12q_Boundry.R` — chr12q boundary **[hpc-first-run]**

```bash
cd "$REPO/analysis/Chr12qEvent"
Rscript --no-save getChr12q_Boundry.R 2>&1 | tail -5
git --no-pager diff --stat -- chr12qEvent.txt data/chr12qEvent.csv 2>/dev/null
git checkout -- chr12qEvent.txt 2>/dev/null
```

PASS: completes and writes `chr12qEvent.txt`.

---

## Step 8 — Sample tables (Table 1 / Supp Tables 1,2) **[verified-local]**

```bash
cd "$REPO/reports"
Rscript --no-save getSampleTable.R 2>&1 | tail -5
ls -l sampleTableProgressionV10.xlsx
Rscript --no-save mkJoinTbl.R 2>&1 | tail -5
ls -l sampleTableProgressionV10__ClinJoin.xlsx
# Both workbooks contain CRDB-derived clinical data and MUST NOT be committed
# (they are gitignored). Remove them after the check so verify leaves no PHI
# artifact behind.
rm -f sampleTableProgressionV10.xlsx sampleTableProgressionV10__ClinJoin.xlsx
```

PASS: both run; `getSampleTable.R` writes a workbook with Summary / AssayTable / Full.CRDB /
Excluded sheets; `mkJoinTbl.R` consumes it and writes the ClinJoin workbook.

---

## Step 9 — Optional / environment-specific

These are not figure-generating and depend on institutional paths. Verify only if you need
to regenerate submission metadata.

### 9a. Submission metadata generators **[hpc-only, known non-portable path]**

- `Geo/makeMasterManifest.R` — GEO GSE244163 manifest.
- `DbGap/Project_04610_WES/mkDbGapMetaData.R` — **reads an absolute path outside the repo**
  (`../../../../CRDB.LIMS/dumps/2022.04.09/...csv`). It only regenerates dbGaP metadata; it
  does not affect any figure/table. Run only if that CRDB dump is reachable from the clone.

---

## Step 10 — Final integrity sweep and cleanup

After the per-step `git checkout --` resets, the working tree should be back to the clean
branch state (no regenerated outputs left over):

```bash
cd "$REPO"
git status -s            # ideally empty (or only intended manifest/README/edits if not yet committed)
```

Confirm no kept script still points at a pruned location (should print nothing):

```bash
grep -rEn 'Pass1/|raw/OldSets/|raw/CGH/rae/(Rae|funcs|getEvent|plot|compute)|GeneVsGeneCorr|doGSEA|makeSampleTable' \
  --include=*.R analysis VennTable figures tables data reports Geo DbGap
```

PASS: empty output.

---

## Sign-off checklist

| # | Check | Result |
|---|---|---|
| 1 | Branch `master`, `~/.Rprofile` helpers present | ☐ |
| 2 | All packages OK incl. `openxlsx` and `bedr`(+bedtools) | ☐ |
| 3a| No broken `data` symlinks | ☐ |
| 3b| All input blobs (`data/`, `data/db/`, `Rlib/`, `FEAT.file`) present | ☐ |
| 3c| Loader smoke test: 349 / 19933×187 / 24904×187 / 20032×24 | ☐ |
| 3d| `data/maf.R` loads (record row counts) | ☐ |
| 4a| `mkVennTable.R` → join table reproduces (text diff trivial) | ☐ |
| 4c| `heatmapV2.R` → Fig 1A | ☐ |
| 5a| `plotRAEProfile.R` → Fig 2A | ☐ |
| 5b| `plotChrRegion.R` → Fig 3B/4A | ☐ |
| 5c| `getCGHEventTableSelected.R` → Supp 4 | ☐ |
| 5d| `getRAEGeneTableSelect.R` → Supp 5 | ☐ |
| 5e| `getRAEGeneTableByBands.R` → Supp 13/14 (needs bedr) | ☐ |
| 6a| WD integration → Supp 7 (dim 127×7); `Rlib/` paths OK | ☐ |
| 6b| DD integration → Supp 12/13 | ☐ |
| 7a| `findProgressionBlocks.R` → Fig 6A inputs | ☐ |
| 7b| `cghJUN.R` → JUN freq (WD ~5.9% / DD ~27%) | ☐ |
| 7c| `getChr12q_Boundry.R` → chr12q event | ☐ |
| 8 | Sample tables (`getSampleTable.R`, `mkJoinTbl.R`) | ☐ |
| 9 | (optional) submission metadata generators | ☐ |
| 10| Final sweep clean; no references to pruned paths | ☐ |

## Not reproducible from this repo (do NOT expect outputs)

Per the manuscript, these were produced **externally** and have no code here (see
`README.md` / `CORE_CODE_MANIFEST.md`):

1. Survival / Kaplan-Meier / Cox — Fig 2B, 4B, Supp Figs 2, 8 (biostatistics).
2. Mutation **Figure 5** + Supp Tables 11/12/15 — the custom-capture validation cohort, via
   cBioPortal (`lipo_cbe_singers_4610_fi`). `data/maf.R` here covers only the 9-patient WES
   **test** set (Supp Tables 9/10).
3. Upstream pipelines — RNA-seq (STAR/HTSeq/DESeq) and WES variant calling
   (`soccin/BIC-variants_pipeline`, `soccin/Variant-PostProcess`); aCGH RAE calling (its
   output is precompiled in `data/db/`).
