peptides <- paste0("LibA_001_", 1:3)
annot <- make_annotation(peptides, pro_id = c("A", "A", "B"))

test_that("annotation columns land between the id and the samples", {
  d <- make_matrix(peptides, 2L)
  out <- annotate_table(d, annot)
  expect_equal(names(out), c("u_pep_id", ANNOT_FIELDS, "sample1", "sample2"))
  expect_equal(out$u_pep_id, peptides)
  expect_equal(out$sample1, d$sample1)
})

test_that("peptide-specific fields are dropped for protein tables", {
  pro <- data.table::data.table(pro_id = c("A", "B"), s1 = c(1, 2))
  pannot <- data.table::data.table(pro_id = c("A", "B"))
  for (f in ANNOT_FIELDS) data.table::set(pannot, j = f, value = c("x", "y"))
  out <- annotate_table(pro, pannot)
  expect_equal(names(out),
               c("pro_id", "UniProt_acc", "taxon_genus", "taxon_species",
                 "gene_symbol", "product", "s1"))
  expect_false(any(c("pep_id", "pos_start", "pos_end", "pep_aa") %in% names(out)))
})

test_that("peptides absent from the annotation get NA, not a shifted row", {
  d <- make_matrix(peptides, 1L)
  partial <- annot[u_pep_id != peptides[2L]]
  out <- annotate_table(d, partial)
  expect_equal(nrow(out), 3L)
  expect_true(is.na(out$gene_symbol[2L]))
  expect_equal(out$gene_symbol[c(1L, 3L)], partial$gene_symbol)
})

test_that("a missing join column is an error", {
  d <- make_matrix(peptides, 1L)
  expect_error(annotate_table(d, annot[, .(pro_id)]), "no 'u_pep_id' column")
})

test_that("a missing annotation field is an error", {
  d <- make_matrix(peptides, 1L)
  expect_error(annotate_table(d, annot[, .(u_pep_id, pro_id)]),
               "missing requested field")
})

test_that("a value missing in the annotation is written as the literal NA", {
  # The annotation has a row for this peptide, but pos_end is empty. That has to
  # come out as text so the column round-trips the same way it always has.
  d <- make_matrix(peptides, 1L)
  a <- data.table::copy(annot)
  data.table::set(a, j = "pos_end", value = c(56L, NA_integer_, 56L))
  out <- annotate_table(d, a)
  expect_identical(out$pos_end, c("56", "NA", "56"))
})

test_that("a peptide with no annotation row keeps a true NA", {
  d <- make_matrix(peptides, 1L)
  partial <- annot[u_pep_id != peptides[2L]]
  out <- annotate_table(d, partial)
  expect_true(is.na(out$pos_end[2L]))
  expect_false(identical(out$pos_end[2L], "NA"))
})

test_that("the two missing cases are distinguishable in one table", {
  d <- make_matrix(peptides, 1L)
  a <- annot[u_pep_id != peptides[3L]]
  data.table::set(a, j = "gene_symbol", value = c("G1", NA_character_))
  out <- annotate_table(d, a)
  expect_identical(out$gene_symbol[2L], "NA")  # present, value missing
  expect_true(is.na(out$gene_symbol[3L]))      # absent from annotation
})
