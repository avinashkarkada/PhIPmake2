write_pair <- function(a, b) {
  da <- tempfile("a"); db <- tempfile("b")
  dir.create(da); dir.create(db)
  data.table::fwrite(a, file.path(da, "x.tsv"), sep = "\t")
  data.table::fwrite(b, file.path(db, "x.tsv"), sep = "\t")
  list(a = da, b = db)
}

base <- data.table::data.table(id = c("p1", "p2"), s1 = c(1, 2), s2 = c(3, 4))

test_that("identical directories compare clean", {
  p <- write_pair(base, base)
  res <- compare_outputs(p$a, p$b, verbose = FALSE)
  expect_equal(res$verdict, "identical")
  expect_equal(res$n_diff_cells, 0L)
})

test_that("differing values are located and quantified", {
  other <- data.table::copy(base)
  other$s2 <- c(3, 9)
  p <- write_pair(base, other)
  res <- compare_outputs(p$a, p$b, verbose = FALSE)
  expect_equal(res$verdict, "values_differ")
  expect_equal(res$n_diff_cells, 1L)
  expect_equal(res$max_abs_diff, 5)
  expect_match(res$detail, "s2\\[p2\\]")
})

test_that("row and column count changes are reported as shape differences", {
  p <- write_pair(base, base[1L])
  expect_equal(compare_outputs(p$a, p$b, verbose = FALSE)$verdict, "shape_differs")

  p2 <- write_pair(base, base[, .(id, s1)])
  expect_equal(compare_outputs(p2$a, p2$b, verbose = FALSE)$verdict, "shape_differs")
})

test_that("reordered identifiers are flagged rather than silently diffed", {
  p <- write_pair(base, base[c(2L, 1L)])
  res <- compare_outputs(p$a, p$b, verbose = FALSE)
  expect_equal(res$verdict, "ids_differ")
})

test_that("files present on only one side are reported", {
  p <- write_pair(base, base)
  data.table::fwrite(base, file.path(p$b, "extra.tsv"), sep = "\t")
  res <- compare_outputs(p$a, p$b, verbose = FALSE)
  expect_true("extra_in_new" %in% res$verdict)
})

test_that("tolerance suppresses floating point noise", {
  other <- data.table::copy(base)
  other$s1 <- base$s1 + 1e-12
  p <- write_pair(base, other)
  expect_equal(compare_outputs(p$a, p$b, verbose = FALSE)$verdict, "identical")
  expect_equal(compare_outputs(p$a, p$b, tolerance = 1e-15, verbose = FALSE)$verdict,
               "values_differ")
})

test_that("byte-identical files are reported as such", {
  p <- write_pair(base, base)
  res <- compare_outputs(p$a, p$b, verbose = FALSE)
  expect_equal(res$bytes, "same")
})

test_that("same values written differently are flagged as formatting only", {
  # Identical data, but one side quotes the header and the other does not.
  da <- tempfile(); db <- tempfile(); dir.create(da); dir.create(db)
  data.table::fwrite(base, file.path(da, "x.tsv"), sep = "\t", quote = TRUE)
  data.table::fwrite(base, file.path(db, "x.tsv"), sep = "\t", quote = FALSE)
  res <- compare_outputs(da, db, verbose = FALSE)
  expect_equal(res$verdict, "identical")   # values agree
  expect_equal(res$bytes, "differ")        # bytes do not
})

test_that("write_matrix quotes the header and character values, not numerics", {
  f <- file.path(tempfile("d"), "out.tsv")
  write_matrix(data.table::data.table(id = "A_001_x", lab = "1", n = 10L), f)
  lines <- readLines(f)
  expect_identical(lines[1], "\"id\"\t\"lab\"\t\"n\"")
  expect_identical(lines[2], "\"A_001_x\"\t\"1\"\t10")
})
