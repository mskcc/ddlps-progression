#!/usr/bin/env bash
#
# verify.sh - automated runner for the protocol in VERIFY.md
# (bash 3.2+ compatible: no associative arrays, so it runs on macOS too)
#
# Confirms the pruned, paper-only public repo (master) is intact and reproduces
# the manuscript's computational outputs. Requires R (>= 4.2.2) with openxlsx
# and bedr, and the bedtools binary on PATH. The legacy ~/.Rprofile helpers
# (cc, len, DATE, suppress, write.xls, write_xlsx) are bundled in R/helpers.R
# and are sourced by every script that needs them via per-directory
# R -> ../R symlinks, so no user dotfile is required.
#
# (Historical: the pruning work happened on branch manu/v_2026, which is now
# the master branch of the public repo.)
#
# USAGE
#   ./verify.sh                 # run every step, print a summary table
#   ./verify.sh 3c 4a 6a        # run only the named steps
#   ./verify.sh --list          # list step ids and what they check
#   ./verify.sh --no-reset      # do NOT git-checkout regenerated outputs (inspect them)
#   ./verify.sh --help
#
# DESIGN
#   - Never aborts on the first failure: every check runs, results are tallied,
#     and a sign-off table is printed at the end. Exit code is the number of
#     hard FAILs (0 = all good).
#   - Working tree safety: any tracked output file a step regenerates is reset
#     with `git checkout --` at the end (unless --no-reset), but ONLY if it was
#     clean before this run. Files you had already modified are left untouched.
#   - Statuses: PASS (clean), FAIL (hard problem - count toward exit code),
#     WARN (needs a human eyeball, e.g. binary xlsx/pdf byte diff, or an
#     optional env dep missing), SKIP (not requested / dependency absent).
#
# This script is read-only with respect to source: it runs scripts and then
# restores their regenerated outputs. It does not edit any .R file.

# ---- not `set -e`: we deliberately tolerate per-check failures ----
set -u
set -o pipefail

# ----------------------------------------------------------------------------
# Locate repo root (script lives at repo root)
# ----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$SCRIPT_DIR}"
cd "$REPO" || { echo "Cannot cd to REPO=$REPO"; exit 99; }

RSCRIPT="${RSCRIPT:-Rscript}"

# Minimum supported R version (floor). Development was on 4.5.1; 4.2.2 is the
# tested floor. Override with: R_MIN=4.3.0 ./verify.sh
R_MIN="${R_MIN:-4.2.2}"

# ----------------------------------------------------------------------------
# Colors (disabled if not a tty)
# ----------------------------------------------------------------------------
if [ -t 1 ]; then
  C_RED=$'\033[31m'; C_GRN=$'\033[32m'; C_YEL=$'\033[33m'
  C_BLU=$'\033[34m'; C_BOLD=$'\033[1m'; C_OFF=$'\033[0m'
else
  C_RED=""; C_GRN=""; C_YEL=""; C_BLU=""; C_BOLD=""; C_OFF=""
fi

# ----------------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------------
DO_RESET=1
REQUESTED=()
for a in "$@"; do
  case "$a" in
    --help|-h)
      sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    --list)
      LIST_ONLY=1 ;;
    --no-reset)
      DO_RESET=0 ;;
    *)
      REQUESTED+=("$a") ;;
  esac
done

# ----------------------------------------------------------------------------
# Result accounting
# ----------------------------------------------------------------------------
declare -a RESULT_ID RESULT_STATUS RESULT_MSG
N_PASS=0; N_FAIL=0; N_WARN=0; N_SKIP=0

record() {  # record <id> <STATUS> <message>
  RESULT_ID+=("$1"); RESULT_STATUS+=("$2"); RESULT_MSG+=("$3")
  case "$2" in
    PASS) N_PASS=$((N_PASS+1)); printf "  %s[PASS]%s %s\n" "$C_GRN" "$C_OFF" "$3" ;;
    FAIL) N_FAIL=$((N_FAIL+1)); printf "  %s[FAIL]%s %s\n" "$C_RED" "$C_OFF" "$3" ;;
    WARN) N_WARN=$((N_WARN+1)); printf "  %s[WARN]%s %s\n" "$C_YEL" "$C_OFF" "$3" ;;
    SKIP) N_SKIP=$((N_SKIP+1)); printf "  %s[SKIP]%s %s\n" "$C_BLU" "$C_OFF" "$3" ;;
  esac
}

step_hdr() { printf "\n%s== Step %s — %s ==%s\n" "$C_BOLD" "$1" "$2" "$C_OFF"; }

# Decide whether a step id was requested (no args = all)
want() {
  [ "${#REQUESTED[@]}" -eq 0 ] && return 0
  local id="$1" r
  for r in ${REQUESTED[@]+"${REQUESTED[@]}"}; do [ "$r" = "$id" ] && return 0; done
  return 1
}

# ----------------------------------------------------------------------------
# Working-tree snapshot: remember which tracked files were ALREADY dirty, so we
# never reset something the user was editing. We only reset files we regenerate
# that were clean at start. (bash 3.2 has no associative arrays -> use a
# newline-delimited string membership test.)
# ----------------------------------------------------------------------------
DIRTY_LIST=""
TO_RESET_LIST=""
snapshot_dirty() {
  DIRTY_LIST="$(git -C "$REPO" diff --name-only; git -C "$REPO" diff --name-only --staged)"
}
was_dirty() {  # was_dirty <path> -> rc 0 if user already had it modified
  printf '%s\n' "$DIRTY_LIST" | grep -Fxq "$1"
}
mark_reset() {  # queue tracked output files (clean-at-start) for end-of-run reset
  local f
  for f in "$@"; do
    f="${f#./}"
    was_dirty "$f" || TO_RESET_LIST="$TO_RESET_LIST$f"$'\n'
  done
}
do_resets() {
  if [ "$DO_RESET" -eq 0 ]; then
    echo; echo "--no-reset: leaving regenerated files in place."; return
  fi
  [ -z "$TO_RESET_LIST" ] && return
  echo; echo "Cleaning up regenerated outputs..."
  printf '%s' "$TO_RESET_LIST" | sort -u | while IFS= read -r f; do
    [ -z "$f" ] && continue
    # Tracked outputs: restore to committed baseline.
    # Untracked outputs (e.g. PHI-bearing sample-table xlsx that must not live
    # in the repo): remove so the verify run leaves no artifact behind.
    if git -C "$REPO" ls-files --error-unmatch -- "$f" >/dev/null 2>&1; then
      git -C "$REPO" checkout -- "$f" 2>/dev/null && echo "  reset    $f"
    else
      rm -f "$REPO/$f" && echo "  removed  $f"
    fi
  done
}

# ----------------------------------------------------------------------------
# Helpers for individual checks
# ----------------------------------------------------------------------------

# Run an R script from its own dir; capture rc and output, and stamp a marker
# file just before launch so we can tell which outputs were FRESHLY written
# (mtime newer than the marker) versus the committed baseline that was already
# on disk.  Usage: run_r <dir> <script>  -> sets RC and ROUT; RUN_MARKER set.
RUN_MARKER=""
trap '[ -n "$RUN_MARKER" ] && rm -f "$RUN_MARKER"' EXIT
run_r() {
  local dir="$1" scr="$2"
  [ -n "$RUN_MARKER" ] && rm -f "$RUN_MARKER"
  RUN_MARKER="$(mktemp "${TMPDIR:-/tmp}/verify_marker.XXXXXX")"
  # ensure the marker's mtime is strictly before the run starts
  sleep 1
  ROUT="$( cd "$REPO/$dir" && "$RSCRIPT" --no-save "$scr" 2>&1 )"
  RC=$?
}

# True if <file> was written/modified after the last run_r marker (a fresh
# write by the script just run), not merely the pre-existing committed file.
wrote_fresh() {  # wrote_fresh <abs_path>
  [ -n "$RUN_MARKER" ] || return 1
  [ -f "$1" ] || return 1
  [ "$1" -nt "$RUN_MARKER" ]
}

# Count of changed content lines in a tracked file's diff (float noise still
# counts here; the human confirms >0 is trivial). 0 = perfect reproduction.
text_diff_changed_lines() {
  git -C "$REPO" --no-pager diff --no-color -- "$1" \
    | grep -E '^[+-]' | grep -Ev '^(\+\+\+|---)' | wc -l | tr -d ' '
}

# Float-tolerant comparison of a tracked text file against its committed
# baseline. Splits each line on whitespace/tabs into cells and compares cell
# by cell: numeric cells must agree to within rel-tol 1e-9; non-numeric cells
# must match exactly. Echoes the number of lines that DIFFER MATERIALLY (i.e.
# changes that are not pure trailing-digit float noise). 0 = reproduces.
# Falls back to raw `text_diff_changed_lines` on any internal error.
text_diff_changed_lines_float_tol() {
  local path="$1" rc out
  out="$("$RSCRIPT" --no-save -e '
args   <- commandArgs(trailingOnly=TRUE)
relpath<- args[1]
reltol <- 1e-9
cur    <- tryCatch(readLines(relpath, warn=FALSE), error=function(e) NULL)
base   <- tryCatch(suppressWarnings(system2("git", c("show", paste0(":",relpath)),
                                            stdout=TRUE, stderr=FALSE)),
                   error=function(e) NULL)
if (is.null(cur) || is.null(base)) { cat("NA"); quit(status=0) }
n      <- max(length(cur), length(base))
length(cur)  <- n
length(base) <- n
cells_equal <- function(a, b) {
  if (identical(a, b)) return(TRUE)
  na <- suppressWarnings(as.numeric(a))
  nb <- suppressWarnings(as.numeric(b))
  if (is.na(na) || is.na(nb)) return(FALSE)
  if (na == 0 && nb == 0) return(TRUE)
  abs(na - nb) <= reltol * max(abs(na), abs(nb))
}
material <- 0L
for (i in seq_len(n)) {
  la <- cur[[i]];  lb <- base[[i]]
  if (is.na(la) && is.na(lb)) next
  if (is.na(la) || is.na(lb)) { material <- material + 1L; next }
  if (identical(la, lb)) next
  ca <- strsplit(la, "\t", fixed=TRUE)[[1]]
  cb <- strsplit(lb, "\t", fixed=TRUE)[[1]]
  if (length(ca) != length(cb)) { material <- material + 1L; next }
  diff_cells <- !mapply(cells_equal, ca, cb)
  if (any(diff_cells)) material <- material + 1L
}
cat(material)
' "$path" 2>/dev/null)"
  rc=$?
  if [ "$rc" -ne 0 ] || [ -z "$out" ] || [ "$out" = "NA" ]; then
    text_diff_changed_lines "$path"
  else
    echo "$out"
  fi
}

# Identify which optional dependency a failed run was missing, by inspecting the
# error text.  Echoes one of: xlsx | bedr | tidygenomics | "" (none detected).
missing_dep() {  # missing_dep <output>
  local o="$1"
  if echo "$o" | grep -qiE "no package called .?xlsx|could not find function .?write\.xlsx"; then
    echo xlsx
  elif echo "$o" | grep -qiE "no package called .?bedr|could not find function .?bedr\."; then
    echo bedr
  elif echo "$o" | grep -qiE "no package called .?tidygenomics|could not find function .?genome_"; then
    echo tidygenomics
  else
    echo ""
  fi
}

# ============================================================================
# STEP DEFINITIONS
# Each is a function step_<id> that runs its checks via record().
# ============================================================================

step_0() {
  step_hdr 0 "Environment (R >= $R_MIN)"
  if ! command -v "$RSCRIPT" >/dev/null 2>&1; then
    record 0 FAIL "Rscript not on PATH (install R or set RSCRIPT=/path/to/Rscript)"
    return
  fi
  # Numeric comparison via R (handles 4.2.2 vs 4.10.x correctly, unlike globbing).
  local rver cmp
  rver="$("$RSCRIPT" -e 'cat(as.character(getRversion()))' 2>/dev/null)"
  cmp="$("$RSCRIPT" -e "cat(as.integer(getRversion() >= '$R_MIN'))" 2>/dev/null)"
  if [ "$cmp" = "1" ]; then
    record 0 PASS "R $rver (>= floor $R_MIN)"
  else
    record 0 FAIL "R $rver is below the floor $R_MIN - upgrade the R module"
  fi
}

step_1() {
  step_hdr 1 "Branch + R/helpers.R"
  local br ref; br="$(git -C "$REPO" branch --show-current)"
  # PASS on master itself, or on any feature branch descended from it (i.e.
  # master's tip is an ancestor of HEAD). Refs are resolved locally; if no
  # local master ref exists, fall back to origin/master.
  if [ "$br" = "master" ]; then
    record 1 PASS "on branch master ($(git -C "$REPO" log --oneline -1))"
  else
    if git -C "$REPO" rev-parse --verify --quiet master >/dev/null; then
      ref="master"
    elif git -C "$REPO" rev-parse --verify --quiet origin/master >/dev/null; then
      ref="origin/master"
    else
      ref=""
    fi
    if [ -n "$ref" ] && git -C "$REPO" merge-base --is-ancestor "$ref" HEAD 2>/dev/null; then
      record 1 PASS "on branch '$br' (descended from $ref - $(git -C "$REPO" log --oneline -1))"
    else
      record 1 WARN "branch is '$br' (not a descendant of master)"
    fi
  fi
  if [ -f "$REPO/R/helpers.R" ]; then
    local ok
    ok="$( cd "$REPO" && "$RSCRIPT" -e 'source("R/helpers.R"); cat(all(sapply(c("cc","write.xls","write_xlsx","DATE","len","suppress"), exists)))' 2>/dev/null)"
    if [ "$ok" = "TRUE" ]; then
      record 1 PASS "R/helpers.R defines cc/write.xls/write_xlsx/DATE/len/suppress"
    else
      record 1 FAIL "R/helpers.R missing one or more helpers (got: '$ok')"
    fi
  else
    record 1 FAIL "R/helpers.R not found"
  fi
}

step_2() {
  step_hdr 2 "Required R packages"
  local out
  # All .xlsx writes in the figure pipeline now use openxlsx (the xlsx
  # migration; see commit 8200088 + reports follow-up). xlsx is therefore
  # no longer in the required list. The CRDB site-info update script
  # (data/raw/CRDB/updateCRDBWithSiteInfo.R) is the last remaining xlsx
  # caller and is NOT exercised by any verify step.
  out="$("$RSCRIPT" -e '
pkgs <- c("tidyverse","limma","edgeR","gplots","IRanges","GenomicRanges",
          "data.table","readxl","openxlsx","digest","RSQLite","stringr",
          "org.Hs.eg.db","AnnotationDbi","knitr",
          "bedr","tidygenomics")
inst <- rownames(installed.packages())
miss <- pkgs[!(pkgs %in% inst)]
cat(if(length(miss)) paste(miss, collapse=" ") else "NONE")
' 2>/dev/null)"
  if [ "$out" = "NONE" ]; then
    record 2 PASS "all required packages installed (openxlsx, bedr, tidygenomics, org.Hs.eg.db, AnnotationDbi, knitr)"
  else
    # bedr/tidygenomics missing = WARN (degrades specific steps); core = FAIL
    local hard=""
    for p in $out; do
      case "$p" in bedr|tidygenomics) ;; *) hard="$hard $p" ;; esac
    done
    if [ -n "$hard" ]; then
      record 2 FAIL "missing core packages:$hard"
    else
      record 2 WARN "missing optional:$out (bedr->bands; tidygenomics->getRAEGeneTableSelect)"
    fi
  fi
  if command -v bedtools >/dev/null 2>&1; then
    record 2 PASS "bedtools on PATH: $(bedtools --version 2>&1 | head -1)"
  else
    record 2 WARN "bedtools not on PATH (bedr step 5e will fail without it)"
  fi
}

step_3a() {
  step_hdr 3a "All 'data' symlinks resolve"
  local broken=0 l
  while IFS= read -r l; do
    if [ ! -e "$l/" ]; then echo "    BROKEN: $l -> $(readlink "$l")"; broken=1; fi
  done < <(find "$REPO" -type l -name data -not -path '*/.git/*')
  if [ "$broken" -eq 0 ]; then
    record 3a PASS "all data symlinks resolve"
  else
    record 3a FAIL "one or more data symlinks broken (see above)"
  fi
}

step_3b() {
  step_hdr 3b "Required input blobs present"
  local missing=0 f
  local files=(
    data/u133a.rda data/cellLines.rda data/targetsUnion.Rdata
    data/sampleTable.csv data/chr12qEvent.csv
    "data/raw/CRDB/FullDataExport_Nick_Socci_04_05_2016___PaperFREEZE_2017_05_03.txt"
    data/raw/CGH/cghCBSRdataFiles
    data/raw/CGH/rae/FEAT.file
    data/db/XX_lumiNorm_20160516_.Rdata
    data/db/CragoProgressionC-lesions.Rdata
    data/db/GeneMatrix.txt
    data/db/CragoProgressionC_EventMatrixV2_20160419_.txt
    data/db/Proj_04610_manu___SOMATIC_FACETS.vep.filtered.maf.gz
    data/db/Proj_3704_Merge_GeneCounts.txt.gz
    data/db/targetsUnion.Rdata
    data/db/maf_colClasses data/db/maf_colnames
    data/db/progressionEventTable.csv
    data/db/averageSignal_3Sort_M3_.txt data/db/human.hg18.genome
    data/db/progressionSet.txt data/db/CNV.file
    analysis/mRNAvsCGH/Rlib/annotation.R analysis/mRNAvsCGH/Rlib/hgu133a.sqlite
  )
  for f in "${files[@]}"; do
    [ -e "$REPO/$f" ] || { echo "    MISSING $f"; missing=$((missing+1)); }
  done
  if [ "$missing" -eq 0 ]; then
    record 3b PASS "all ${#files[@]} input blobs present (incl. promoted Rlib/ + FEAT.file)"
  else
    record 3b FAIL "$missing input blob(s) missing (see above)"
  fi
}

step_3c() {
  step_hdr 3c "Loader smoke test (expected dims)"
  local out
  out="$( cd "$REPO/data" && "$RSCRIPT" --no-save -e '
source("R/helpers.R")
suppressWarnings(suppressMessages({
  source("sampleTable.R"); source("progressionSet.R"); source("crdb.R")
  source("cghUBM.R");      source("cghGeneMatrix.R");   source("rnaSeq.R")
  source("miRNATargets.R");source("shRNA.R");           source("hg18.R")
  source("regions.R")
}))
cat("ST",nrow(sampleTable),"UBM",paste(dim(ubm$a0),collapse="x"),
    "GCGH",paste(dim(geneCGH),collapse="x"),"RNA",paste(dim(rnaSeq$ds),collapse="x"),"\n")
' 2>&1 )"
  if echo "$out" | grep -q "ST 349 UBM 19933x187 GCGH 24904x187 RNA 20032x24"; then
    record 3c PASS "loaders OK: 349 / 19933x187 / 24904x187 / 20032x24"
  else
    record 3c FAIL "loader dims differ or error: $(echo "$out" | tail -2 | tr '\n' ' ')"
  fi
}

step_3d() {
  step_hdr 3d "data/maf.R loads WES MAF"
  local out
  out="$( cd "$REPO/data" && "$RSCRIPT" --no-save -e '
source("R/helpers.R")
suppressWarnings(suppressMessages(source("maf.R")))
cat("COMPLETE",nrow(mafs$complete),"PATIENT",nrow(mafs$patient),"\n")
' 2>&1 )"
  if echo "$out" | grep -qE "COMPLETE [0-9]+ PATIENT [0-9]+"; then
    record 3d PASS "maf loads: $(echo "$out" | grep -oE 'COMPLETE [0-9]+ PATIENT [0-9]+')"
  else
    record 3d FAIL "maf load failed: $(echo "$out" | tail -2 | tr '\n' ' ')"
  fi
}

# Generic "run script + check a TEXT output reproduces" check.
# Judges the TEXT file even when the run later dies at an xlsx write, because
# the .txt is produced first.
check_text_repro() {  # <id> <dir> <script> <committed_txt>
  local id="$1" dir="$2" scr="$3" txt="$4"
  run_r "$dir" "$scr"
  mark_reset "$dir/$txt"
  if ! wrote_fresh "$REPO/$dir/$txt"; then
    local dep; dep="$(missing_dep "$ROUT")"
    if [ -n "$dep" ]; then
      record "$id" WARN "$scr needs '$dep' - $txt not (re)written; install the missing R package"
    else
      record "$id" FAIL "$scr did not write $txt (rc=$RC): $(echo "$ROUT" | tail -1)"
    fi
    return
  fi
  local changed; changed="$(text_diff_changed_lines "$dir/$txt")"
  if [ "$changed" = "0" ]; then
    record "$id" PASS "$scr -> $txt reproduces (no diff)"
    return
  fi
  # Raw diff found something; rule out pure float noise (cross-platform
  # IEEE-754 trailing-digit rounding) before flagging as WARN.
  local material; material="$(text_diff_changed_lines_float_tol "$dir/$txt")"
  if [ "$material" = "0" ]; then
    record "$id" PASS "$scr -> $txt reproduces ($changed line(s) differ only in trailing-digit float noise, content identical)"
  else
    record "$id" WARN "$scr -> $txt has $material materially changed line(s) (of $changed raw) - inspect"
  fi
}

# Generic "run script + check BINARY output(s) got FRESHLY written" check.
# Distinguishes a real fresh write from the committed baseline already on disk.
check_bin_written() {  # <id> <dir> <script> <committed_bin...>
  local id="$1" dir="$2" scr="$3"; shift 3
  run_r "$dir" "$scr"
  local b nfresh=0 ntotal=0
  for b in "$@"; do
    mark_reset "$dir/$b"
    ntotal=$((ntotal+1))
    wrote_fresh "$REPO/$dir/$b" && nfresh=$((nfresh+1))
  done
  if [ "$nfresh" -eq "$ntotal" ] && [ "$RC" -eq 0 ]; then
    record "$id" PASS "$scr completed; $nfresh/$ntotal output(s) freshly written"
  elif [ "$nfresh" -gt 0 ] && [ "$RC" -ne 0 ]; then
    record "$id" WARN "$scr rc=$RC; $nfresh/$ntotal outputs written - inspect: $(echo "$ROUT" | tail -1)"
  else
    local dep; dep="$(missing_dep "$ROUT")"
    if [ -n "$dep" ]; then
      record "$id" WARN "$scr needs '$dep' - output not written; install the missing R package"
    else
      record "$id" FAIL "$scr produced no fresh output (rc=$RC): $(echo "$ROUT" | tail -1)"
    fi
  fi
}

step_4a() {
  step_hdr 4a "Venn v14 (Fig 1B + Supp 3)"
  check_text_repro 4a VennTable mkVennTable.R joinTableCragoProgression_v14_.txt
  mark_reset VennTable/cghCellLine_miRNA_Venn_v14_.pdf \
             VennTable/mrnaGeneSetsVenn_v14_.pdf \
             VennTable/vennDiagrams_mRNA_v14_.pdf \
             VennTable/joinTableCragoProgression_v14_.xlsx
}

step_4c() {
  step_hdr 4c "heatmapV2.R (Fig 1A)"
  if [ ! -e "$REPO/figures/mRNAHeatmap/joinTableCragoProgression_v14_.txt" ]; then
    record 4c FAIL "join-table symlink in figures/mRNAHeatmap does not resolve"
    return
  fi
  check_bin_written 4c figures/mRNAHeatmap heatmapV2.R heatMap_v14.pdf
}

step_5a() {
  step_hdr 5a "plotRAEProfile.R (Fig 2A / Supp 1)"
  run_r figures/CGHProfiles plotRAEProfile.R
  local f nfresh=0
  for f in "$REPO"/figures/CGHProfiles/raePlotA0D0_*.pdf "$REPO"/figures/CGHProfiles/raePlotA0D0_*.png; do
    [ -e "$f" ] || continue
    mark_reset "${f#$REPO/}"
    wrote_fresh "$f" && nfresh=$((nfresh+1))
  done
  if [ "$RC" -eq 0 ] && [ "$nfresh" -ge 1 ]; then
    record 5a PASS "plotRAEProfile.R freshly wrote $nfresh raePlotA0D0_* file(s)"
  else
    local dep; dep="$(missing_dep "$ROUT")"
    if [ -n "$dep" ]; then
      record 5a WARN "plotRAEProfile.R needs '$dep' - not written"
    else
      record 5a FAIL "plotRAEProfile.R rc=$RC, $nfresh fresh file(s): $(echo "$ROUT" | tail -1)"
    fi
  fi
}

step_5b() {
  step_hdr 5b "plotChrRegion.R (Fig 3B / 4A)"
  check_bin_written 5b figures/CGHProfiles plotChrRegion.R \
    "chr6Region_WD+DD_v4.pdf" chr6RegionTest_v3.pdf chr6RegionTest_WDOnly_v3.pdf
}

step_5c() {
  step_hdr 5c "getCGHEventTableSelected.R (Supp 4)"
  check_bin_written 5c tables/CGHEventTable getCGHEventTableSelected.R \
    cghEventTable___Manuscript__SelectRegions1.xlsx
}

step_5d() {
  step_hdr 5d "getRAEGeneTableSelect.R (Supp 5)"
  check_bin_written 5d tables/CGHEventTable getRAEGeneTableSelect.R \
    cghGeneTable___Manuscript___SelectRegions1.xlsx
}

step_5e() {
  step_hdr 5e "getRAEGeneTableByBands.R (Supp 13/14, needs bedr)"
  check_bin_written 5e tables/CGHEventTable getRAEGeneTableByBands.R \
    raeGeneTableSigRegionsV1.xlsx
}

step_6a() {
  step_hdr 6a "WD integration (Fig 3A / Supp 7) + Rlib repoint"
  if ! grep -q 'source("Rlib/annotation.R")' "$REPO/analysis/mRNAvsCGH/integrate_cgh_mrna__WDonly__v2.R"; then
    record 6a FAIL "WDonly script does not source Rlib/annotation.R (repoint missing?)"
    return
  fi
  if grep -q 'Pass1/' "$REPO/analysis/mRNAvsCGH/integrate_cgh_mrna__WDonly__v2.R"; then
    record 6a FAIL "WDonly script still references Pass1/ path"
    return
  fi
  check_bin_written 6a analysis/mRNAvsCGH integrate_cgh_mrna__WDonly__v2.R \
    suppTable_4_with_U133A__WDpEvent_vs_WD_without_v2.xlsx \
    diff_mRNA_vs_CopyNumber_by_blocks___WDonly_v7_.rda
  # stronger: dim check 127x7 (only if file present and readxl available)
  if [ -f "$REPO/analysis/mRNAvsCGH/suppTable_4_with_U133A__WDpEvent_vs_WD_without_v2.xlsx" ]; then
    local d; d="$("$RSCRIPT" -e 'suppressMessages(a<-readxl::read_xlsx("analysis/mRNAvsCGH/suppTable_4_with_U133A__WDpEvent_vs_WD_without_v2.xlsx")); cat(nrow(a),ncol(a))' 2>/dev/null)"
    [ "$d" = "127 7" ] && record 6a PASS "WD suppTable_4 dim 127x7 (matches reference)" \
                        || record 6a WARN "WD suppTable_4 dim '$d' (ref 127 7) - inspect"
  fi
}

step_6b() {
  step_hdr 6b "DD integration (Fig 6B / Supp 8,13) + Rlib repoint"
  if ! grep -q 'source("Rlib/annotation.R")' "$REPO/analysis/mRNAvsCGH/integrate_cgh_mrna__DDonly__v2.R"; then
    record 6b FAIL "DDonly script does not source Rlib/annotation.R"
    return
  fi
  if grep -q 'Pass1/' "$REPO/analysis/mRNAvsCGH/integrate_cgh_mrna__DDonly__v2.R"; then
    record 6b FAIL "DDonly script still references Pass1/ path"
    return
  fi
  check_bin_written 6b analysis/mRNAvsCGH integrate_cgh_mrna__DDonly__v2.R \
    suppTable_12_with_U133A__DDpEvent_vs_DD_without_v3.xlsx \
    diff_mRNA_vs_CopyNumber_by_blocks___DDonly_v8_.rda
}

step_7a() {
  step_hdr 7a "findProgressionBlocks.R (Fig 6A inputs)"
  # writes TWO workbooks via openxlsx: __Manifest.xlsx (early) and the main
  # _v1_.xlsx (later). Track both so neither leaks into the tree.
  check_bin_written 7a analysis/ProgressionBlocks findProgressionBlocks.R \
    progressionBlocks_WDvsDD_Expr_v1_.xlsx \
    progressionBlocks_WDvsDD_Expr_v1__Manifest.xlsx
}

step_7b() {
  step_hdr 7b "cghJUN.R (JUN 1p32 frequency)"
  run_r analysis/JUN_Freq cghJUN.R
  if [ "$RC" -eq 0 ]; then
    record 7b PASS "cghJUN.R ran (expect WDLS gain 5/88 ~5.7% / DDLS 25/93 ~27%, p<0.001)"
    # show from the Supp Table 14 headline through the full breakdown
    echo "$ROUT" | sed -n '/Supplemental Table 14/,$p' | sed 's/^/      | /'
  else
    record 7b FAIL "cghJUN.R rc=$RC: $(echo "$ROUT" | tail -1)"
  fi
}

step_7c() {
  step_hdr 7c "getChr12q_Boundry.R (chr12q boundary)"
  run_r analysis/Chr12qEvent getChr12q_Boundry.R
  mark_reset analysis/Chr12qEvent/chr12qEvent.txt
  if [ "$RC" -eq 0 ] && [ -f "$REPO/analysis/Chr12qEvent/chr12qEvent.txt" ]; then
    record 7c PASS "getChr12q_Boundry.R -> chr12qEvent.txt written"
  else
    record 7c FAIL "getChr12q_Boundry.R rc=$RC or no output: $(echo "$ROUT" | tail -1)"
  fi
}

step_8() {
  step_hdr 8 "Sample tables (Table 1 / Supp 1,2)"
  check_bin_written 8a reports getSampleTable.R sampleTableProgressionV10.xlsx
  # mkJoinTbl consumes the workbook getSampleTable just wrote
  check_bin_written 8b reports mkJoinTbl.R sampleTableProgressionV10__ClinJoin.xlsx
}

step_10() {
  step_hdr 10 "Final sweep: no kept script points at pruned paths"
  local hits
  hits="$(grep -rEn 'Pass1/|raw/OldSets/|raw/CGH/rae/(Rae|funcs|getEvent|plot|compute)|GeneVsGeneCorr|doGSEA|makeSampleTable' \
        --include=*.R analysis VennTable figures tables data reports 2>/dev/null)"
  if [ -z "$hits" ]; then
    record 10 PASS "no references to pruned locations in kept scripts"
  else
    record 10 FAIL "kept scripts reference pruned paths:"
    echo "$hits" | sed 's/^/      /'
  fi
}

# ----------------------------------------------------------------------------
# Step registry (ordered)
# ----------------------------------------------------------------------------
ALL_STEPS=(0 1 2 3a 3b 3c 3d 4a 4c 5a 5b 5c 5d 5e 6a 6b 7a 7b 7c 8 10)

# Description lookup (bash 3.2-safe: case statement, not associative array)
step_desc() {
  case "$1" in
    0)  echo "R present (>= floor)" ;;
    1)  echo "branch + R/helpers.R" ;;
    2)  echo "required R packages (openxlsx, bedr, tidygenomics, org.Hs.eg.db, AnnotationDbi, knitr)" ;;
    3a) echo "data symlinks resolve" ;;
    3b) echo "input blobs present" ;;
    3c) echo "loader smoke test (dims)" ;;
    3d) echo "maf.R loads" ;;
    4a) echo "Venn v14 -> join table reproduces" ;;
    4c) echo "heatmapV2.R -> Fig 1A" ;;
    5a) echo "plotRAEProfile.R -> Fig 2A" ;;
    5b) echo "plotChrRegion.R -> Fig 3B/4A" ;;
    5c) echo "getCGHEventTableSelected.R -> Supp 4" ;;
    5d) echo "getRAEGeneTableSelect.R -> Supp 5" ;;
    5e) echo "getRAEGeneTableByBands.R -> Supp 13/14 (bedr)" ;;
    6a) echo "WD integration -> Supp 7 (127x7) + Rlib" ;;
    6b) echo "DD integration -> Supp 8/13 + Rlib" ;;
    7a) echo "findProgressionBlocks.R -> Fig 6A inputs" ;;
    7b) echo "cghJUN.R -> JUN freq" ;;
    7c) echo "getChr12q_Boundry.R -> chr12q" ;;
    8)  echo "sample tables (getSampleTable/mkJoinTbl)" ;;
    8a) echo "getSampleTable.R -> Table 1 / Supp 1,2" ;;
    8b) echo "mkJoinTbl.R -> clinical join table" ;;
    10) echo "final sweep: no pruned-path refs" ;;
    *)  echo "$1" ;;
  esac
}

if [ "${LIST_ONLY:-0}" = "1" ]; then
  echo "Step ids and checks:"
  for s in "${ALL_STEPS[@]}"; do printf "  %-4s %s\n" "$s" "$(step_desc "$s")"; done
  exit 0
fi

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
printf "%sverify.sh%s  REPO=%s\n" "$C_BOLD" "$C_OFF" "$REPO"
printf "Branch: %s   %s\n" "$(git -C "$REPO" branch --show-current)" "$(date)"
if [ "${#REQUESTED[@]}" -gt 0 ]; then
  printf "Running only: %s\n" "${REQUESTED[*]}"
fi

snapshot_dirty

for s in "${ALL_STEPS[@]}"; do
  if want "$s"; then
    "step_${s//./_}" 2>&1 || record "$s" FAIL "step crashed unexpectedly"
  fi
done

do_resets

# ----------------------------------------------------------------------------
# Summary table (one row per step id, collapsed to its WORST sub-result so a
# step with two sub-checks does not print twice). Severity FAIL>WARN>SKIP>PASS.
# ----------------------------------------------------------------------------
sev() { case "$1" in FAIL) echo 3;; WARN) echo 2;; SKIP) echo 1;; *) echo 0;; esac; }

printf "\n%s================ SUMMARY ================%s\n" "$C_BOLD" "$C_OFF"
printf "%-5s %-7s %s\n" "STEP" "RESULT" "CHECK"
printf -- "-------------------------------------------------\n"
n_results="${#RESULT_ID[@]}"
# distinct ids in first-seen order
seen_ids=""
i=0
while [ "$i" -lt "$n_results" ]; do
  the_id="${RESULT_ID[$i]}"
  case " $seen_ids " in *" $the_id "*) ;; *) seen_ids="$seen_ids $the_id" ;; esac
  i=$((i+1))
done
for the_id in $seen_ids; do
  # worst status across all records for this id
  worst="PASS"; j=0
  while [ "$j" -lt "$n_results" ]; do
    if [ "${RESULT_ID[$j]}" = "$the_id" ]; then
      [ "$(sev "${RESULT_STATUS[$j]}")" -gt "$(sev "$worst")" ] && worst="${RESULT_STATUS[$j]}"
    fi
    j=$((j+1))
  done
  case "$worst" in
    PASS) col="$C_GRN" ;; FAIL) col="$C_RED" ;; WARN) col="$C_YEL" ;; *) col="$C_BLU" ;;
  esac
  printf "%-5s %s%-7s%s %s\n" "$the_id" "$col" "$worst" "$C_OFF" "$(step_desc "$the_id")"
done
printf -- "-------------------------------------------------\n"
printf "%sPASS %d   %sFAIL %d   %sWARN %d   %sSKIP %d%s\n" \
  "$C_GRN" "$N_PASS" "$C_RED" "$N_FAIL" "$C_YEL" "$N_WARN" "$C_BLU" "$N_SKIP" "$C_OFF"
echo
echo "WARN = needs a human eyeball (binary byte-diff, float noise, or optional dep missing)."
echo "See VERIFY.md for the PASS/FAIL criteria behind each step."

exit "$N_FAIL"
