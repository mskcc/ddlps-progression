# How the Paired Log Fold Change Is Computed

## The simple approach: mean of within-pair differences

The most intuitive way to compute a paired log fold change would be:

1. For each of the 4 patient pairs, compute the difference in log expression between their WD and NF sample
2. Average those 4 differences

```
Pair 1:  log(WD₁) − log(NF₁) = d₁
Pair 2:  log(WD₂) − log(NF₂) = d₂
Pair 3:  log(WD₃) − log(NF₃) = d₃
Pair 4:  log(WD₄) − log(NF₄) = d₄

Log fold change = (d₁ + d₂ + d₃ + d₄) / 4
```

This works and is easy to understand. The GLM approach does something very similar, but with one important refinement.

## What the GLM does differently: weighted averaging

The core difference is that **the GLM doesn't treat all four pairs as equally informative**. Instead it weights each pair by how reliably it can measure expression.

Here is the intuition: RNA-seq counts are noisy in a specific way — genes expressed at very low levels are measured imprecisely (a count of 3 vs 4 is a huge relative difference), while highly expressed genes are measured much more reliably (a count of 3000 vs 3001 is negligible noise). The negative binomial model captures this formally.

So the GLM computes something closer to:

```
Log fold change = weighted average of (d₁, d₂, d₃, d₄)

where pairs with higher, more reliable counts get more weight
```

## Why it matters for a paired design

The `~Patient+Tissue` design matrix tells the model: *"first account for the fact that patients differ from each other overall, then estimate the WD-vs-NF difference."* The Patient terms soak up all the between-patient variability (some patients may just have higher expression overall), so the Tissue estimate is based purely on **within-pair differences** — exactly as you'd want in a paired analysis.

If the four pairs were equally precise, the GLM result and the simple mean of within-pair differences would be **numerically identical**. The GLM is preferred because it handles the unequal precision of count data in a statistically principled way, and it naturally handles the paired structure without requiring a separate preprocessing step.

## One-sentence summary

The log fold change was estimated as the `TissueWD` coefficient from a negative binomial GLM with patient as a blocking factor, equivalent to a precision-weighted average of within-pair log expression differences across the four matched WD/NF pairs.
