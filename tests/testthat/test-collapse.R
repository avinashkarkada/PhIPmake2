peptides <- paste0("LibA_001_", 1:6)
data <- data.table::data.table(u_pep_id = peptides, s1 = c(10, 20, 30, 40, 50, 60))
# p3 is deliberately absent from the annotation.
annot <- data.table::data.table(
  u_pep_id = peptides[-3L],
  pro_id   = c("A", "A", "B", "B", "B"))

test_that("map_peptides_to_proteins preserves row alignment", {
  pro <- suppressWarnings(map_peptides_to_proteins(data$u_pep_id, annot))
  expect_length(pro, nrow(data))
  expect_equal(pro, c("A", "A", NA, "B", "B", "B"))
})

test_that("unmapped peptides are reported", {
  expect_warning(map_peptides_to_proteins(data$u_pep_id, annot),
                 "1 of 6 peptides have no annotation")
})

test_that("prosum groups the correct peptides when an annotation is missing", {
  # Dropping the NA from the protein map first would shorten it to length 5,
  # which R recycles against the 6 data rows and yields A = 90, B = 120.
  pro <- suppressWarnings(map_peptides_to_proteins(data$u_pep_id, annot))
  out <- collapse_peptides(data, pro, "sum")

  expect_equal(out$pro_id, c("A", "B"))
  expect_equal(out$s1, c(30, 150))
  expect_false(identical(out$s1, c(90, 120)))
})

test_that("promax takes the maximum within each protein", {
  pro <- suppressWarnings(map_peptides_to_proteins(data$u_pep_id, annot))
  out <- collapse_peptides(data, pro, "max")
  expect_equal(out$pro_id, c("A", "B"))
  expect_equal(out$s1, c(20, 60))
})

test_that("a mis-sized protein vector is rejected instead of recycled", {
  short <- stats::na.omit(suppressWarnings(
    map_peptides_to_proteins(data$u_pep_id, annot)))
  expect_error(collapse_peptides(data, as.character(short), "sum"),
               "Pass the full-length vector")
})

test_that("protein order follows first appearance among annotated peptides", {
  d <- data.table::data.table(u_pep_id = paste0("LibA_001_", 1:4), s1 = 1:4)
  a <- data.table::data.table(u_pep_id = d$u_pep_id, pro_id = c("Z", "Y", "Z", "Y"))
  out <- collapse_peptides(d, map_peptides_to_proteins(d$u_pep_id, a), "sum")
  expect_equal(out$pro_id, c("Z", "Y"))
  expect_equal(out$s1, c(4, 6))
})

test_that("multiple sample columns are collapsed independently", {
  d <- data.table::data.table(u_pep_id = paste0("LibA_001_", 1:4),
                              s1 = c(1, 2, 3, 4), s2 = c(40, 30, 20, 10))
  a <- data.table::data.table(u_pep_id = d$u_pep_id, pro_id = c("P", "P", "Q", "Q"))
  pro <- map_peptides_to_proteins(d$u_pep_id, a)
  expect_equal(collapse_peptides(d, pro, "sum")$s2, c(70, 30))
  expect_equal(collapse_peptides(d, pro, "max")$s1, c(2, 4))
})

test_that("na_rows='drop' discards peptides with any missing value", {
  d <- data.table::data.table(u_pep_id = paste0("LibA_001_", 1:3),
                              s1 = c(1, NA, 3), s2 = c(10, 20, 30))
  a <- data.table::data.table(u_pep_id = d$u_pep_id, pro_id = c("P", "P", "P"))
  pro <- map_peptides_to_proteins(d$u_pep_id, a)
  # Peptide 2 has an NA in s1, so the whole row is discarded from the protein.
  expect_equal(collapse_peptides(d, pro, "sum", na_rows = "drop")$s2, 40)
  # Keeping the row uses every available value per sample instead.
  expect_equal(collapse_peptides(d, pro, "sum", na_rows = "keep")$s2, 60)
})

test_that("collapsing with no annotated peptides yields an empty table", {
  pro <- rep(NA_character_, nrow(data))
  out <- collapse_peptides(data, pro, "sum")
  expect_equal(nrow(out), 0L)
  expect_equal(names(out), c("pro_id", "s1"))
})

test_that("the input matrix is not modified in place", {
  before <- data.table::copy(data)
  pro <- suppressWarnings(map_peptides_to_proteins(data$u_pep_id, annot))
  collapse_peptides(data, pro, "sum")
  expect_equal(data, before)
})
