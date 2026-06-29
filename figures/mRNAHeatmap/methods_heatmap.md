# Figure 1A heatmap — methods, terminology, and accuracy

This documents what `heatmapV2.R` does, why the figure looks the way it does, and
how the reported classification numbers are defined — so this never has to be
re-derived. (It cost real time once; read this instead.)

## What the script does

`heatmapV2.R`:
1. Loads U133A expression (`data(u133a)`) and the sample table (`data(sampleTable)`).
2. Selects all U133A samples of TYPE `NF` (normal fat) or `WD` (WDLS) →
   **90 entities** (62 WD, 28 NF).
3. Uses the gene set in `joinTableCragoProgression_v14_.txt` (the DE signature from
   the v14 Venn analysis) as the features — column `PROBE.U133A`.
4. Row-wise z-scores the expression matrix (center + scale per gene), clips to ±4.
5. Draws `gplots::heatmap.2(t(dZZc), dendro="row", ...)` → clusters the **rows**,
   which are the **samples**, and color-bars them by true class (NF = blue, WD = red).
6. Cuts that same dendrogram at **k=3** and writes the accuracy numbers to
   `heatmap_accuracy_v14.txt`.

Outputs: `heatMap_v14.pdf` (2 pages: same heatmap, sample-ID labels then class
labels) and `heatmap_accuracy_v14.txt`.

## Why k=3 (not k=2)

The published Figure 1A reads **three** brackets off the dendrogram:

| Figure bracket | Composition | Call |
|---|---|---|
| "21/22 normal fat" (top) | 22 entities: 21 NF, 1 WD | NF |
| "35/42 WDLS" (middle) | 42 entities: 7 NF, 35 WD | WD |
| "26/26 WDLS" (bottom) | 26 entities: 0 NF, 26 WD | WD |

Total 22+42+26 = 90. These three brackets are **exactly reproduced** by the script's
clustering (Euclidean distance + complete linkage, the `heatmap.2` defaults) cut at
**k=3**. Do **not** cut at k=2: the top split of this tree is unbalanced (64 vs 26)
and both top branches are majority-WD, so k=2 collapses to a meaningless 68.9%. The
NF samples form their own pure sub-branch that only appears at k≥3 — which is why the
figure shows three groups, not two.

## Clustering parameters (locked)

- distance: **Euclidean** (`dist` default) on the per-gene z-scored, ±4-clipped matrix
- linkage: **complete** (`hclust` default)
- These are `heatmap.2`'s defaults and they reproduce the published brackets exactly.
  Other linkages (ward.D2, average) give cleaner top-level splits and ~96.7% accuracy
  but do **not** match the published figure's 21/22 · 35/42 · 26/26 brackets. The
  figure was made with the defaults; keep them.

## Terminology — do NOT call this "supervised" or "unsupervised"

This is the trap that starts reviewer holy wars:

- The **clustering algorithm** is unsupervised (hierarchical clustering never sees the
  labels). But the **features** (the genes) were chosen *using* the labels — they are
  the WD-vs-NF differentially expressed genes. So:
  - Calling it **"unsupervised"** invites the selection-bias critique (Ambroise &
    McLachlan 2002; Simon et al. 2003): "you picked genes that already know the
    answer, of course it separates."
  - Calling it **"supervised clustering"** invites the purist objection that
    clustering is unsupervised by definition (the phrase reads as an oxymoron).

**Resolution:** describe the two steps separately and avoid both adjectives:

> "Hierarchical clustering of samples using the genes differentially expressed
>  between WDLS and normal fat..."

This states the supervision honestly (the feature set is the DE signature) while
naming the actual, unimpeachable algorithm (hierarchical clustering).

## Accuracy / sensitivity / specificity / precision (the numbers)

Each of the three brackets is labeled by its majority class → predicted class.
Positive class = **WDLS (tumor)**; negative = normal fat. Confusion matrix:

```
                truth WD   truth NF
predicted WD       61          7
predicted NF        1         21
```
TP=61, FP=7, FN=1, TN=21, N=90.

| Metric | Definition | Plain English | Value |
|---|---|---|---|
| Accuracy | (TP+TN)/N | fraction of all 90 called correctly | **82/90 = 91.1%** |
| Sensitivity (recall, TPR) | TP/(TP+FN) | of true WDLS, fraction caught | **61/62 = 98.4%** |
| Specificity (TNR) | TN/(TN+FP) | of true normal fat, fraction caught | **21/28 = 75.0%** |
| Precision (PPV) | TP/(TP+FP) | of those called WDLS, fraction correct | **61/68 = 89.7%** |

### Where the old "98% accuracy" came from
The manuscript previously said "98% accuracy." That number is actually the
**sensitivity** (61/62 = 98.4%): only one WDLS sample landed in the normal-fat
bracket. It was the right number with the wrong word. The honest, defensible
statement gives all three:

> "...grouped WDLS and normal fat with 91% accuracy (82/90 samples; 98% sensitivity
>  and 75% specificity for WDLS)."

The signature is excellent at *catching* tumors (98% sensitivity) but weaker at
*ruling them out* (75% specificity — 7 of 28 normal-fat samples drift into a WDLS
bracket; visible as red stripes in the upper region of the heatmap).

## Manuscript wording to keep in sync

Two places in `Crago_TAB2_with_Figures__NDS-Claude.docx` (supplement needs none):
1. Results, "Dysregulation of gene expression in WDLS" paragraph — the
   "supervised clustering ... 98% accuracy" sentence.
2. Figure 1 caption, panel A — "supervised clustering".

Both should use "hierarchical clustering of samples using the DE gene set" and the
91% / 98% / 75% numbers above. Regenerate the numbers any time with:

```
cd figures/mRNAHeatmap && Rscript --no-save heatmapV2.R
cat heatmap_accuracy_v14.txt
```
