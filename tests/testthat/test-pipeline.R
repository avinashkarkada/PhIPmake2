# End-to-end run on a miniature two-library screen, asserting the output
# contract: file names, per-library directories, and pan files that are the
# concatenation of their sub-library parts.

build_screen <- function(n_samples = 3L) {
  wd <- file.path(tempfile("screen"))
  md <- file.path(tempfile("metadata"))
  dir.create(wd, recursive = TRUE); dir.create(md, recursive = TRUE)

  libs <- c("LibAlpha_001", "LibBeta_001")
  peps <- unlist(lapply(libs, make_peptides, n = 6L))
  pro <- c(rep(c("A1", "A2"), each = 3L), rep(c("B1", "B2", "B3"), each = 2L))

  set.seed(42)
  n <- length(peps)
  counts <- make_matrix(peps, n_samples,
                        values = matrix(as.numeric(sample(1:100, n * n_samples, TRUE)),
                                        nrow = n))
  fc <- make_matrix(peps, n_samples,
                    values = matrix(round(runif(n * n_samples, 0.5, 8), 3), nrow = n))
  enr <- make_matrix(peps, n_samples,
                     values = matrix(round(runif(n * n_samples), 4), nrow = n))
  hitm <- make_matrix(peps, n_samples,
                      values = matrix(sample(0:1, n * n_samples, TRUE), nrow = n))

  data.table::fwrite(counts, file.path(wd, "counts.csv"))
  data.table::fwrite(fc, file.path(wd, "fold_change.csv"))
  data.table::fwrite(enr, file.path(wd, "enrichment.csv"))
  data.table::fwrite(hitm, file.path(wd, "Hits.csv"))

  for (l in libs) {
    keep <- startsWith(peps, paste0(l, "_"))
    write_metadata(md, l, make_annotation(peps[keep], pro[keep]),
                   pairs = data.table::data.table(V1 = peps[keep][1L],
                                                  V2 = peps[keep][2L]))
  }
  list(wd = wd, md = md, libs = libs, peptides = peps, n_samples = n_samples)
}

test_that("a full run produces the expected file set", {
  s <- build_screen()
  expect_silent(suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = c("counts", "foldchange", "enrichment", "hits"),
                  screen_name = "screenX", metadata_path = s$md, verbose = FALSE))))

  expected_pan <- c(
    "screenX_Counts.tsv", "screenX_Counts_annotated.tsv",
    "screenX_Counts_prosum.tsv", "screenX_Counts_prosum_annotated.tsv",
    "screenX_FoldChange.tsv", "screenX_FoldChange_annotated.tsv",
    "screenX_EdgeR.tsv", "screenX_EdgeR_annotated.tsv",
    "screenX_EdgeR_promax.tsv", "screenX_EdgeR_promax_annotated.tsv",
    "screenX_Hits.tsv", "screenX_Hits_annotated.tsv",
    "screenX_Hits_prosum.tsv", "screenX_Hits_prosum_annotated.tsv",
    "screenX_Hits_counts.tsv", "screenX_Hits_counts_annotated.tsv",
    "screenX_Hits_foldchange.tsv", "screenX_Hits_foldchange_annotated.tsv",
    "screenX_Hits_foldchange_promax.tsv", "screenX_Hits_foldchange_promax_annotated.tsv",
    "screenX_Hits_enrichment.tsv", "screenX_Hits_enrichment_annotated.tsv",
    "screenX_Hits_enrichment_promax.tsv", "screenX_Hits_enrichment_promax_annotated.tsv")
  for (f in expected_pan)
    expect_true(file.exists(file.path(s$wd, f)), info = f)

  for (l in s$libs) {
    expect_true(dir.exists(file.path(s$wd, l)))
    expect_true(file.exists(file.path(s$wd, l, sprintf("screenX_%s_Counts.tsv", l))))
    expect_true(file.exists(file.path(s$wd, l, sprintf("screenX_%s_EdgeR_promax_annotated.tsv", l))))
  }
})

test_that("pan files are exactly the concatenation of their sub-library files", {
  s <- build_screen()
  suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = c("counts", "hits"), screen_name = "screenX",
                  metadata_path = s$md, verbose = FALSE)))

  for (stem in c("Counts_prosum", "Hits_prosum")) {
    pan <- data.table::fread(file.path(s$wd, sprintf("screenX_%s.tsv", stem)))
    subs <- data.table::rbindlist(lapply(s$libs, function(l)
      data.table::fread(file.path(s$wd, l, sprintf("screenX_%s_%s.tsv", l, stem)))))
    expect_equal(pan, subs, info = stem)
  }
})

test_that("sub-library files partition the peptides without loss or overlap", {
  s <- build_screen()
  suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = "counts", screen_name = "screenX",
                  metadata_path = s$md, verbose = FALSE)))

  ids <- unlist(lapply(s$libs, function(l)
    data.table::fread(file.path(s$wd, l, sprintf("screenX_%s_Counts.tsv", l)))[[1L]]))
  expect_setequal(ids, s$peptides)
  expect_equal(anyDuplicated(ids), 0L)
})

test_that("hit masking is applied with the right default per data type", {
  s <- build_screen()
  suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = c("counts", "foldchange", "enrichment", "hits"),
                  screen_name = "screenX", metadata_path = s$md, verbose = FALSE)))

  hits <- data.table::fread(file.path(s$wd, "Hits.csv"))
  hfc <- data.table::fread(file.path(s$wd, "screenX_Hits_foldchange.tsv"))
  hen <- data.table::fread(file.path(s$wd, "screenX_Hits_enrichment.tsv"))
  data.table::setkeyv(hits, names(hits)[1L]); data.table::setkeyv(hfc, names(hfc)[1L])
  data.table::setkeyv(hen, names(hen)[1L])

  zero <- hits[["sample1"]] == 0
  expect_true(all(hfc[["sample1"]][zero] == 1))   # fold change non-hits -> 1
  expect_true(all(hen[["sample1"]][zero] == 0))   # enrichment non-hits -> 0
})

test_that("polyclonal output is produced and bounded by the hit count", {
  s <- build_screen()
  suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = c("hits", "polyclonal"), screen_name = "screenX",
                  metadata_path = s$md, verbose = FALSE)))

  poly <- data.table::fread(file.path(s$wd, "screenX_Polyclonal.tsv"))
  prosum <- data.table::fread(file.path(s$wd, "screenX_Hits_prosum.tsv"))
  expect_equal(poly$pro_id, prosum$pro_id)
  # An epitope count can never exceed the number of hit peptides.
  expect_true(all(poly$sample1 <= prosum$sample1))
})

test_that("the dash-separated stage string from the original pipeline works", {
  s <- build_screen()
  suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = "Counts-Hits", screen_name = "screenX",
                  metadata_path = s$md, verbose = FALSE)))
  expect_true(file.exists(file.path(s$wd, "screenX_Counts_prosum.tsv")))
  expect_true(file.exists(file.path(s$wd, "screenX_Hits_prosum.tsv")))
  expect_false(file.exists(file.path(s$wd, "screenX_EdgeR.tsv")))
})

test_that("an unknown stage is rejected", {
  s <- build_screen()
  expect_error(run_phipmake2(s$wd, stages = "Bogus", metadata_path = s$md),
               "Unknown stage")
})

test_that("unimplemented stages are skipped with a warning, not an error", {
  s <- build_screen()
  expect_warning(
    suppressMessages(run_phipmake2(
      s$wd, stages = "Counts-AVARDA-Epitopefindr", screen_name = "screenX",
      metadata_path = s$md, verbose = FALSE)),
    "not implemented and skipped: avarda, epitopefindr")
  # The runnable part of the plan still ran.
  expect_true(file.exists(file.path(s$wd, "screenX_Counts_prosum.tsv")))
})

test_that("a plan of only unimplemented stages is an error", {
  s <- build_screen()
  expect_error(
    suppressWarnings(run_phipmake2(s$wd, stages = "AVARDA", metadata_path = s$md)),
    "No runnable stages")
})

test_that("drake_params.tsv is still honoured", {
  s <- build_screen()
  data.table::fwrite(data.table::data.table(
    params = c("screen_name", "enrichment_type", "metadata_path"),
    value  = c("fromParams", "EdgeR", s$md)),
    file.path(s$wd, "drake_params.tsv"), sep = "\t")

  suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = "counts", verbose = FALSE)))
  expect_true(file.exists(file.path(s$wd, "fromParams_Counts.tsv")))
})

test_that("a run is idempotent", {
  s <- build_screen()
  run <- function() suppressMessages(suppressWarnings(
    run_phipmake2(s$wd, stages = "counts", screen_name = "screenX",
                  metadata_path = s$md, verbose = FALSE)))
  run()
  first <- readLines(file.path(s$wd, "screenX_Counts_prosum.tsv"))
  run()
  expect_equal(readLines(file.path(s$wd, "screenX_Counts_prosum.tsv")), first)
})
