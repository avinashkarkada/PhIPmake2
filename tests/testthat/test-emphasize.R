peptides <- paste0("LibA_001_", 1:4)
data <- data.table::data.table(u_pep_id = peptides,
                               s1 = c(2.5, 3.5, 4.5, 5.5),
                               s2 = c(9, 8, 7, 6))
hits <- data.table::data.table(u_pep_id = peptides,
                               s1 = c(1, 0, 1, 0),
                               s2 = c(0, 0, 1, 1))

test_that("non-hit cells are replaced with the default", {
  out <- emphasize_hits(data, hits, 0)
  expect_equal(out$s1, c(2.5, 0, 4.5, 0))
  expect_equal(out$s2, c(0, 0, 7, 6))
})

test_that("the fold-change default of 1 is honoured", {
  out <- emphasize_hits(data, hits, 1)
  expect_equal(out$s1, c(2.5, 1, 4.5, 1))
})

test_that("identifiers are preserved and the input is untouched", {
  before <- data.table::copy(data)
  out <- emphasize_hits(data, hits, 0)
  expect_equal(out$u_pep_id, peptides)
  expect_equal(data, before)
})

test_that("mismatched shapes are rejected", {
  expect_error(emphasize_hits(data, hits[1:2], 0), "4 rows but hits has 2")
  expect_error(emphasize_hits(data, hits[, .(u_pep_id, s1)], 0),
               "3 columns but hits has 2")
})

test_that("masking by position is refused when row order differs", {
  shuffled <- hits[c(2L, 1L, 3L, 4L)]
  expect_error(emphasize_hits(data, shuffled, 0), "not in the same row order")
})
