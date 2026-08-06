# Metadata directories are not laid out uniformly. Most annotation files carry a
# date prefix, some do not, directories also hold Specialty and u_pep_id files
# and a Historic/ subdirectory, and the top level contains loose scripts and a
# directory literally named NA.

fake_metadata <- function(entries, extra_dirs = character()) {
  root <- tempfile("md"); dir.create(root)
  for (lib in names(entries)) {
    d <- file.path(root, lib); dir.create(d)
    for (f in entries[[lib]])
      data.table::fwrite(data.table::data.table(u_pep_id = "L_001_A", pro_id = "P1",
                                                gene_symbol = "G"),
                         file.path(d, f), sep = "\t")
    dir.create(file.path(d, "Historic"), showWarnings = FALSE)
  }
  for (d in extra_dirs) dir.create(file.path(root, d), showWarnings = FALSE)
  root
}

test_that("a date-prefixed annotation file is found", {
  md <- fake_metadata(list(LibAlpha_001 = c(
    "20190530_LibAlpha_001_Universal.tsv",
    "20190530_LibAlpha_001_Specialty.tsv",
    "LibAlpha_001_u_pep_id.txt")))
  a <- read_annotation("LibAlpha_001", md)
  expect_true(all(c("u_pep_id", "pro_id") %in% names(a)))
})

test_that("an annotation file with no date prefix is found", {
  md <- fake_metadata(list(LibTheta_001 = "LibTheta_001_Universal.tsv"))
  expect_s3_class(read_annotation("LibTheta_001", md), "data.table")
})

test_that("Specialty files are never mistaken for Universal ones", {
  md <- fake_metadata(list(LibA_001 = c("20200101_LibA_001_Specialty.tsv",
                                        "20200101_LibA_001_Universal.tsv")))
  expect_silent(read_annotation("LibA_001", md))
})

test_that("a merged _000 library resolves to its own directory", {
  md <- fake_metadata(list(LibAlpha_000 = "20200101_LibAlpha_000_Universal.tsv",
                           LibAlpha_001 = "20190530_LibAlpha_001_Universal.tsv"))
  expect_silent(read_annotation("LibAlpha_000", md))
})

test_that("the Historic subdirectory is not treated as a candidate file", {
  md <- fake_metadata(list(LibA_001 = "20200101_LibA_001_Universal.tsv"))
  dir.create(file.path(md, "LibA_001", "Historic", "sub"), recursive = TRUE)
  expect_silent(read_annotation("LibA_001", md))
})

test_that("a library with no metadata directory is an error", {
  md <- fake_metadata(list(LibA_001 = "20200101_LibA_001_Universal.tsv"),
                      extra_dirs = c("NA", "Universal", "Historic"))
  expect_error(read_annotation("LibZ_009", md), "No metadata directory")
})

test_that("a directory with no Universal file is an error", {
  md <- fake_metadata(list(LibA_001 = "20200101_LibA_001_Specialty.tsv"))
  expect_error(read_annotation("LibA_001", md), "No '_Universal' annotation file")
})

test_that("competing Universal files pick the exact library version", {
  md <- fake_metadata(list(LibA_001 = c("20190101_LibA_001_Universal.tsv",
                                        "20200101_LibA_002_Universal.tsv")))
  expect_silent(read_annotation("LibA_001", md))
})

test_that("trimmedpairs files are found alongside the annotation", {
  md <- fake_metadata(list(LibA_001 = c("20200101_LibA_001_Universal.tsv",
                                        "LibA_001_trimmedpairs.tsv")))
  expect_equal(nrow(read_pairs("LibA_001", md)), 1L)
})

test_that("a library with no pairs file returns an empty table", {
  md <- fake_metadata(list(LibA_001 = "20200101_LibA_001_Universal.tsv"))
  p <- read_pairs("LibA_001", md)
  expect_equal(nrow(p), 0L)
  expect_gte(ncol(p), 2L)
})
