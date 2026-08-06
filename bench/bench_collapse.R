#!/usr/bin/env Rscript
# Protein collapse: loop over proteins vs one grouped pass.
#
# The loop version below is the implementation this package replaces, kept here
# verbatim so the two can be run on identical inputs. Checks that they agree
# when the annotation is complete, how they differ when it is not, and how the
# runtimes compare.
#
#   Rscript bench/bench_collapse.R [--full]
#
# --full adds one pass at the scale of a large sub-library (106,679 peptides
# x 145 samples).

suppressPackageStartupMessages({
  library(data.table); library(magrittr); library(phipmake2)
})

full <- "--full" %in% commandArgs(TRUE)

# --------------------------------------------------------------- loop version
promax_loop <- function(data, margin = 2)
  apply(data, margin, function(x) x %>% stats::na.omit() %>% max())
prosum_loop <- function(data, margin = 2)
  apply(data, 2, function(x) x %>% stats::na.omit() %>% sum())

# Unchanged apart from taking the annotation join as an argument, so that both
# versions see identical inputs.
collapse_by_loop <- function(data, annot, fun = c("max", "sum")) {
  fun <- match.arg(fun)
  pro_ids <- annot$pro_id[match(data[, 1], annot$u_pep_id)] %>% stats::na.omit()
  proteins <- unique(pro_ids)

  lib <- data.frame(matrix(nrow = length(proteins), ncol = ncol(data)))
  colnames(lib) <- c("pro_id", colnames(data)[-1])
  lib[, 1] <- proteins

  collapse <- if (fun == "max") promax_loop else prosum_loop
  for (j in seq_along(proteins)) {
    lib[j, -1] <- data[pro_ids == proteins[j], -1] %>% stats::na.omit() %>% collapse()
  }
  lib
}

# ------------------------------------------------------------------- fixtures
make_data <- function(n_pep, n_samp, n_pro, coverage = 1) {
  set.seed(1)
  peptides <- sprintf("LibAlpha_001_%08d", seq_len(n_pep))
  d <- data.frame(u_pep_id = peptides,
                  matrix(round(runif(n_pep * n_samp) * 100, 2), nrow = n_pep),
                  stringsAsFactors = FALSE)
  colnames(d)[-1] <- paste0("sample", seq_len(n_samp))
  keep <- if (coverage >= 1) seq_len(n_pep) else
    sort(sample.int(n_pep, floor(n_pep * coverage)))
  annot <- data.frame(u_pep_id = peptides[keep],
                      pro_id = sprintf("PRO%06d", ((keep - 1L) %% n_pro) + 1L),
                      stringsAsFactors = FALSE)
  list(data = d, annot = annot)
}

timed <- function(expr) {
  invisible(gc(reset = TRUE, full = TRUE))
  t <- system.time(value <- force(expr))[["elapsed"]]
  peak <- sum(gc()[, "max used"] * c(56, 8)) / 1024^2   # MB, cells + vector words
  list(value = value, seconds = t, peak_mb = peak)
}

# ----------------------------------------------------- 1. correctness, complete
cat("\n=== 1. Agreement when every peptide is annotated ===\n")
fx <- make_data(4000L, 12L, 800L, coverage = 1)
loop_out <- collapse_by_loop(fx$data, fx$annot, "max")
pro <- map_peptides_to_proteins(fx$data$u_pep_id, fx$annot)
grouped_out <- collapse_peptides(fx$data, pro, "max")
cat(sprintf("  protein order identical : %s\n", identical(loop_out$pro_id, grouped_out$pro_id)))
cat(sprintf("  values identical        : %s\n",
            isTRUE(all.equal(as.matrix(loop_out[, -1]), as.matrix(grouped_out[, -1]),
                             check.attributes = FALSE))))

# --------------------------------------------------- 2. divergence, incomplete
cat("\n=== 2. Behaviour when some peptides lack annotation ===\n")
fx2 <- make_data(4000L, 4L, 800L, coverage = 0.98)
loop_inc <- collapse_by_loop(fx2$data, fx2$annot, "sum")
pro_inc <- suppressWarnings(map_peptides_to_proteins(fx2$data$u_pep_id, fx2$annot))
grouped_inc <- collapse_peptides(fx2$data, pro_inc, "sum")

# Ground truth computed independently, without either implementation.
truth <- vapply(grouped_inc$pro_id, function(p) {
  rows <- which(!is.na(pro_inc) & pro_inc == p)
  sum(fx2$data$sample1[rows])
}, numeric(1L))

common <- intersect(loop_inc$pro_id, grouped_inc$pro_id)
i1 <- match(common, loop_inc$pro_id); i2 <- match(common, grouped_inc$pro_id)
cat(sprintf("  proteins reported, loop    : %d\n", nrow(loop_inc)))
cat(sprintf("  proteins reported, grouped : %d\n", nrow(grouped_inc)))
cat(sprintf("  loop disagrees with truth on    %d / %d proteins\n",
            sum(abs(loop_inc$sample1[i1] - truth[i2]) > 1e-9), length(common)))
cat(sprintf("  grouped disagrees with truth on %d / %d proteins\n",
            sum(abs(grouped_inc$sample1[i2] - truth[i2]) > 1e-9), length(common)))

# ------------------------------------------------------------- 3. scaling
cat("\n=== 3. Scaling (145 samples, 5 peptides per protein) ===\n")
cat(sprintf("%9s %12s %12s %10s %11s %11s\n",
            "peptides", "loop (s)", "grouped (s)", "speedup", "loop MB", "grouped MB"))
for (n in c(2000L, 4000L, 8000L, 16000L, 32000L)) {
  fx <- make_data(n, 145L, n %/% 5L, coverage = 1)
  a <- timed(collapse_by_loop(fx$data, fx$annot, "max"))
  p <- map_peptides_to_proteins(fx$data$u_pep_id, fx$annot)
  b <- timed(collapse_peptides(fx$data, p, "max"))
  stopifnot(isTRUE(all.equal(as.matrix(a$value[, -1]), as.matrix(b$value[, -1]),
                             check.attributes = FALSE)))
  cat(sprintf("%9d %12.2f %12.3f %9.0fx %11.0f %11.0f\n",
              n, a$seconds, b$seconds, a$seconds / b$seconds, a$peak_mb, b$peak_mb))
  rm(fx, a, b, p); invisible(gc(FALSE))
}

# ------------------------------------------------------------- 4. real scale
if (full) {
  cat("\n=== 4. Large sub-library scale: 106,679 peptides x 145 samples ===\n")
  fx <- make_data(106679L, 145L, 9431L, coverage = 1)
  p <- map_peptides_to_proteins(fx$data$u_pep_id, fx$annot)
  b <- timed(collapse_peptides(fx$data, p, "max"))
  cat(sprintf("  grouped : %.2f s, peak %.0f MB\n", b$seconds, b$peak_mb))
  a <- timed(collapse_by_loop(fx$data, fx$annot, "max"))
  cat(sprintf("  loop    : %.2f s, peak %.0f MB\n", a$seconds, a$peak_mb))
  cat(sprintf("  speedup : %.0fx\n", a$seconds / b$seconds))
  stopifnot(isTRUE(all.equal(as.matrix(a$value[, -1]), as.matrix(b$value[, -1]),
                             check.attributes = FALSE)))
  cat("  results identical: TRUE\n")
}

cat("\n")
