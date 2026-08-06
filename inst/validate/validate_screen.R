#!/usr/bin/env Rscript
# Validate phipmake2 against an existing phipmake output directory.
#
#   Rscript validate_screen.R --screen <dir> --out <dir> [options]
#
#   --screen    Directory holding the input matrices AND the existing phipmake
#               output to compare against. Not written to.
#   --out       Fresh directory for phipmake2 output. Created if absent.
#   --metadata  Peptide library metadata directory.
#               Defaults to $PHIPMAKE2_METADATA
#   --stages    Default Counts-FoldChange-Enrichment-Hits
#   --threads   data.table threads. Default 4.
#   --preflight-only  Run the metadata checks and stop.
#
# Runs three phases: a preflight over the metadata, the run itself with timing
# and peak RSS, then a file-by-file comparison against the existing output.

suppressPackageStartupMessages({library(data.table); library(phipmake2)})

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
flag <- function(f) f %in% args

screen   <- opt("--screen")
outdir   <- opt("--out")
metadata <- opt("--metadata", Sys.getenv("PHIPMAKE2_METADATA", unset = ""))
stages   <- opt("--stages", "Counts-FoldChange-Enrichment-Hits")
threads  <- as.integer(opt("--threads", "4"))

if (is.null(screen) || is.null(outdir))
  stop("--screen and --out are both required.", call. = FALSE)
if (!dir.exists(screen)) stop("No such screen directory: ", screen, call. = FALSE)

rule <- function(x) cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")
peak_rss_mb <- function() {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  l <- grep("^VmHWM:", readLines(f), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", l)) / 1024
}

# ---------------------------------------------------------------- preflight
rule("PHASE 1: preflight")

inputs <- c(counts = "counts.csv", foldchange = "fold_change.csv",
            enrichment = "enrichment.csv", hits = "Hits.csv")
present <- file.exists(file.path(screen, inputs))
for (i in seq_along(inputs))
  cat(sprintf("  %-12s %-18s %s\n", names(inputs)[i], inputs[i],
              if (present[i]) "found" else "MISSING"))
if (!any(present)) stop("No input matrices in ", screen, call. = FALSE)

ref <- file.path(screen, inputs[present][1L])
cat("\n  Reading identifiers from", basename(ref), "...\n")
ids <- fread(ref, select = 1L, showProgress = FALSE)[[1L]]
cat(sprintf("  %s peptides\n", format(length(ids), big.mark = ",")))

libs <- tryCatch(names(peptide_library_index(ids)), error = function(e) {
  cat("\n  FAILED to parse peptide identifiers:\n    ", conditionMessage(e), "\n")
  quit(status = 1L)
})
cat(sprintf("  %d libraries: %s\n\n", length(libs), paste(libs, collapse = ", ")))

need <- c("pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa",
          "taxon_genus", "taxon_species", "gene_symbol", "product")
report <- rbindlist(lapply(libs, function(l) {
  a <- tryCatch(suppressWarnings(read_annotation(l, metadata)),
                error = function(e) conditionMessage(e))
  if (is.character(a))
    return(data.table(library = l, annotation = "ERROR", peptides = NA_integer_,
                      missing_fields = substr(a, 1, 60), unmapped = NA_integer_,
                      pairs = NA_integer_))
  p <- tryCatch(read_pairs(l, metadata), error = function(e) NULL)
  n_lib <- sum(library_of(ids) == l)
  mapped <- suppressWarnings(map_peptides_to_proteins(
    ids[library_of(ids) == l], a, warn_unmapped = FALSE))
  data.table(library = l, annotation = "ok", peptides = n_lib,
             missing_fields = paste(setdiff(c(need, "pro_id"), names(a)), collapse = ","),
             unmapped = sum(is.na(mapped)),
             pairs = if (is.null(p)) NA_integer_ else nrow(p))
}))
print(report)

if (any(report$annotation == "ERROR")) {
  cat("\n  Preflight FAILED: annotation could not be read for one or more libraries.\n")
  quit(status = 1L)
}
if (any(nzchar(report$missing_fields))) {
  cat("\n  Preflight FAILED: annotation is missing required fields.\n")
  quit(status = 1L)
}
cat("\n  Preflight passed.\n")
cat(sprintf("  Peptides with no annotation: %s of %s (%.2f%%)\n",
            format(sum(report$unmapped), big.mark = ","),
            format(sum(report$peptides), big.mark = ","),
            100 * sum(report$unmapped) / sum(report$peptides)))
cat("  These are excluded from protein-level output, and are where results will\n")
cat("  differ from the previous run.\n")

if (flag("--preflight-only")) quit(status = 0L)

# ---------------------------------------------------------------------- run
rule("PHASE 2: run")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (f in inputs[present]) {
  dest <- file.path(outdir, f)
  if (!file.exists(dest)) file.symlink(normalizePath(file.path(screen, f)), dest)
}

t0 <- Sys.time()
run <- tryCatch({
  run_phipmake2(wd = outdir, stages = stages, screen_name = basename(screen),
                metadata_path = metadata, threads = threads, verbose = TRUE)
  TRUE
}, error = function(e) { cat("\n  RUN FAILED: ", conditionMessage(e), "\n"); FALSE })
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (!isTRUE(run)) quit(status = 1L)
cat(sprintf("\n  Elapsed:  %.1f min\n", elapsed))
cat(sprintf("  Peak RSS: %.0f MB\n", peak_rss_mb()))

# ------------------------------------------------------------------ compare
rule("PHASE 3: compare against existing output")

cmp <- compare_outputs(old = screen, new = outdir, verbose = FALSE)
if (!nrow(cmp)) {
  cat("  No comparable files found.\n"); quit(status = 0L)
}

# A processed screen directory holds more than phipmake output: the run
# configuration, and results from later pipeline steps such as ARscore, ARscape
# and VARscore. None of those are ours to reproduce, so set them aside rather
# than reporting them as missing.
foreign <- grepl("drake_params|ARscore|ARScape|VARscore|validation_report",
                 cmp$file, ignore.case = TRUE)
if (any(foreign)) {
  cat("  Ignoring", sum(foreign), "file(s) not produced by phipmake:",
      paste(basename(cmp$file[foreign]), collapse = ", "), "\n")
  cmp <- cmp[!foreign]
}
if (!nrow(cmp)) {
  cat("  Nothing left to compare.\n"); quit(status = 0L)
}

cat("\n  Verdicts:\n")
for (v in names(sort(table(cmp$verdict), decreasing = TRUE)))
  cat(sprintf("    %-16s %d\n", v, sum(cmp$verdict == v)))

cat(sprintf("\n  Byte-identical files: %d of %d\n",
            sum(cmp$bytes == "same", na.rm = TRUE), nrow(cmp)))
fmt_only <- cmp[verdict == "identical" & bytes == "differ"]
if (nrow(fmt_only)) {
  cat(sprintf("  Same values, different formatting: %d\n", nrow(fmt_only)))
  for (i in seq_len(min(nrow(fmt_only), 5L))) cat("    ", fmt_only$file[i], "\n", sep = "")
}

protein <- grepl("prosum|promax|Polyclonal", cmp$file)
cat("\n  Peptide-level files (expected identical):\n")
pep <- cmp[!protein]
for (v in names(sort(table(pep$verdict), decreasing = TRUE)))
  cat(sprintf("    %-16s %d\n", v, sum(pep$verdict == v)))

cat("\n  Protein-level files (may differ where peptides were unannotated):\n")
pro <- cmp[protein]
for (v in names(sort(table(pro$verdict), decreasing = TRUE)))
  cat(sprintf("    %-16s %d\n", v, sum(pro$verdict == v)))

bad <- cmp[verdict != "identical"]
if (nrow(bad)) {
  cat("\n  Details (first 20):\n")
  for (i in seq_len(min(nrow(bad), 20L)))
    cat(sprintf("    %-70s %s\n      %s\n", bad$file[i], bad$verdict[i],
                substr(bad$detail[i], 1, 160)))
}

fwrite(cmp, file.path(outdir, "validation_report.tsv"), sep = "\t")
cat("\n  Full report: ", file.path(outdir, "validation_report.tsv"), "\n", sep = "")

# Only a genuine disagreement counts as a failure. A file present on one side
# and not the other is a difference in what was run, not in what was computed.
unexpected <- pep[verdict %chin% c("values_differ", "shape_differs", "ids_differ")]
if (nrow(unexpected)) {
  cat("\n  ATTENTION: ", nrow(unexpected),
      " peptide-level file(s) differ. These were expected to be identical.\n", sep = "")
  for (i in seq_len(min(nrow(unexpected), 10L)))
    cat("    ", unexpected$file[i], "\n", sep = "")
  quit(status = 2L)
}
cat("\n  All peptide-level files identical.\n")
quit(status = 0L)
