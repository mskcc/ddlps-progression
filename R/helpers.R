#
# Convenience helpers previously provided by the author's personal
# ~/.Rprofile. They are sourced explicitly by each script that uses them
# so the published code does not depend on a user dotfile.
#
# Each script directory that needs these helpers has an `R -> ../R`
# (or `../../R`) symlink alongside its `data` symlink, so every script
# loads this file uniformly with:
#
#   source("R/helpers.R")
#

# paste(..., sep = "_") with a friendlier default.
cc <- function(..., sep = "_", collapse = NULL) {
    paste(..., sep = sep, collapse = collapse)
}

# Short alias for length().
len <- function(x) length(x)

# Compact YYYYMMDD date stamp used in some output filenames.
DATE <- function() format(Sys.Date(), "%Y%m%d")

# Quiet package attach.
suppress <- function(expr) suppressPackageStartupMessages(expr)

# Tab-delimited writer matching the legacy ~/.Rprofile write.xls().
write.xls <- function(dd, filename, row.names = TRUE, col.names = NA, ...) {
    if (!is.data.frame(dd)) {
        dd <- data.frame(dd, check.names = FALSE)
    }
    if (!row.names) {
        col.names <- TRUE
    }
    utils::write.table(dd, file = filename, sep = "\t", quote = FALSE,
                       col.names = col.names, row.names = row.names, ...)
}

# Single-sheet xlsx writer used by a handful of scripts. Wraps openxlsx so
# the package list in README.md stays unchanged.
write_xlsx <- function(x, file, ...) {
    openxlsx::write.xlsx(x, file = file, ...)
}
