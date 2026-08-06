#!/usr/bin/env Rscript
# Command-line entry point for phipmake2.
#
# Usage:
#   Rscript phipmake2.R --wd <screen dir> [options]
#
# Options:
#   --wd PATH             Screen directory containing the input matrices (required)
#   --stages LIST         Dash- or comma-separated stages. Default:
#                         Counts-FoldChange-Enrichment-Hits
#                         Add -Polyclonal to include polyclonal scoring.
#   --screen-name NAME    Output prefix. Default: basename(--wd)
#   --metadata PATH       Peptide library metadata directory
#   --threads N           data.table threads. Default: all available
#   --quiet               Suppress progress output
#
# Exit status is 0 on success, 1 on error.

suppressPackageStartupMessages(library(phipmake2))

args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Option ", flag, " requires a value.", call. = FALSE)
  args[i + 1L]
}
has_flag <- function(flag) flag %in% args

if (has_flag("--help") || has_flag("-h") || length(args) == 0L) {
  cat(readLines(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1L]))[3:20],
      sep = "\n")
  quit(status = 0L)
}

wd <- get_opt("--wd")
if (is.null(wd)) stop("--wd is required.", call. = FALSE)

stages <- get_opt("--stages", "Counts-FoldChange-Enrichment-Hits")
stages <- unlist(strsplit(stages, "[-,]"))

threads <- get_opt("--threads")
if (!is.null(threads)) threads <- as.integer(threads)

ok <- tryCatch({
  run_phipmake2(
    wd            = wd,
    stages        = stages,
    screen_name   = get_opt("--screen-name"),
    metadata_path = get_opt("--metadata"),
    threads       = threads,
    verbose       = !has_flag("--quiet")
  )
  TRUE
}, error = function(e) {
  message("phipmake2 failed: ", conditionMessage(e))
  FALSE
})

quit(status = if (isTRUE(ok)) 0L else 1L)
