test_that("parse_u_pep_id splits base, version and sequence", {
  p <- parse_u_pep_id(c("LibAlpha_001_ATGC", "LibBeta_003_GGGG"))
  expect_equal(p$base, c("LibAlpha", "LibBeta"))
  expect_equal(p$version, c("001", "003"))
  expect_equal(p$sequence, c("ATGC", "GGGG"))
  expect_equal(p$library, c("LibAlpha_001", "LibBeta_003"))
})

test_that("base names containing underscores still parse", {
  # Splitting on "_" and assuming three fields breaks here.
  p <- parse_u_pep_id("My_Lib_002_ATGC")
  expect_equal(p$base, "My_Lib")
  expect_equal(p$version, "002")
  expect_equal(p$sequence, "ATGC")
})

test_that("sequences containing underscores still parse", {
  p <- parse_u_pep_id("LibA_001_ATG_CGT")
  expect_equal(p$base, "LibA")
  expect_equal(p$sequence, "ATG_CGT")
})

test_that("malformed identifiers fail loudly rather than silently", {
  expect_error(parse_u_pep_id("NoVersionHere"), "do not match")
  expect_warning(parse_u_pep_id("NoVersionHere", strict = FALSE), "do not match")
})

test_that("library_of merges multiple versions of one base to _000", {
  ids <- c("LibAlpha_001_A", "LibAlpha_002_B", "LibGamma_001_C")
  expect_equal(library_of(ids),
               c("LibAlpha_000", "LibAlpha_000", "LibGamma_001"))
})

test_that("library_of keeps a single-version base at its own version", {
  expect_equal(library_of(c("LibBeta_003_A", "LibBeta_003_B")),
               c("LibBeta_003", "LibBeta_003"))
})

test_that("library names are not matched by substring", {
  # A substring match would merge these two libraries into one output file.
  ids <- c("LibBeta_001_A", "LibBetaX_001_B")
  expect_equal(library_of(ids), c("LibBeta_001", "LibBetaX_001"))
})

test_that("peptide_library_index returns row indices per library in order", {
  ids <- c("LibB_001_A", "LibA_001_B", "LibB_001_C")
  idx <- peptide_library_index(ids)
  expect_equal(names(idx), c("LibB_001", "LibA_001"))
  expect_equal(idx[["LibB_001"]], c(1L, 3L))
  expect_equal(idx[["LibA_001"]], 2L)
})
