#
# JUN (1p32) aCGH copy-number frequency, WDLS vs DDLS.
#
# Reproduces the "1p32 (JUN) gain" row of Supplemental Table 14 and the value
# cited in the Results ("identified in 27% of DDLS (n=25) and 5.9% of WDLS
# samples (n=5)"):
#
#   Locus           WDLS            DDLS            P value
#   1p32 (JUN) gain 5 (5.9%)        25 (27%)        <0.001
#

source("R/helpers.R")

suppress(library(tidyverse))

# The data() loaders attach Bioconductor packages (IRanges, S4Vectors, ...);
# silence their attach/masking chatter so the report below is readable.
suppressMessages(suppressWarnings({
    data(cghGeneMatrix)
    data(sampleTable)
}))

LOCUS <- "1p32 (JUN) gain"

# Gene-level aCGH call for JUN, restricted to aCGH-profiled samples.
# geneCGH values: <0 loss, 0 neutral, >0 gain.
i_jun <- str_which(rownames(geneCGH), ":JUN$")
bio_id <- sampleTable |> filter(aCGH == "Y") |> pull(BIO_ID)

jun <- tibble(
        BIO_ID = bio_id,
        TYPE   = sampleTable$TYPE[match(bio_id, sampleTable$BIO_ID)],
        call   = as.numeric(geneCGH[i_jun, bio_id])
    ) |>
    mutate(
        TYPE  = factor(TYPE, levels = c("WD", "DD")),
        Event = factor(
                    case_when(call < 0 ~ "Loss", call > 0 ~ "Gain", TRUE ~ "Neutral"),
                    levels = c("Loss", "Neutral", "Gain")
                )
    )

# ---- Full copy-number breakdown (auditable detail) -------------------------

n_by_type <- jun |> count(TYPE, name = "n_total")

detail <- jun |>
    count(TYPE, Event, name = "n", .drop = FALSE) |>
    left_join(n_by_type, by = "TYPE") |>
    mutate(pct = n / n_total) |>
    select(TYPE, Event, n, n_total, pct) |>
    arrange(TYPE, Event)

# ---- Headline: the Supplemental Table 14 row -------------------------------

# Fisher exact test, gain vs. not-gain across WDLS vs. DDLS.
gain_tab <- jun |>
    mutate(Gain = if_else(Event == "Gain", "gain", "not")) |>
    count(Gain, TYPE) |>
    pivot_wider(names_from = TYPE, values_from = n, values_fill = 0) |>
    column_to_rownames("Gain") |>
    as.matrix()

p_value <- fisher.test(gain_tab)$p.value

supp_row <- detail |>
    filter(Event == "Gain") |>
    transmute(
        Locus   = LOCUS,
        Subtype = str_glue("{TYPE}LS"),
        Summary = str_glue("{n} ({round(100 * pct, 1)}%) of {n_total}")
    ) |>
    pivot_wider(names_from = Subtype, values_from = Summary) |>
    mutate(P_value = if_else(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)))

# ---- Report ---------------------------------------------------------------

cat("Supplemental Table 14 row -- JUN (1p32) copy-number gain, WDLS vs DDLS\n")
cat("(Fisher exact test, gain vs. not-gain):\n\n")
print(as.data.frame(supp_row), row.names = FALSE)

cat("\nFull copy-number breakdown:\n\n")
detail |>
    mutate(Frequency = str_glue("{n}/{n_total} ({round(100 * pct, 1)}%)")) |>
    select(Subtype = TYPE, Event, Frequency) |>
    as.data.frame() |>
    print(row.names = FALSE)
