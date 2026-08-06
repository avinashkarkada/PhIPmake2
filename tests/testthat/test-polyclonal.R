peptides <- paste0("LibA_001_", 1:4)
annot <- data.table::data.table(u_pep_id = peptides,
                                pro_id = c("P", "P", "P", "Q"))
pro <- map_peptides_to_proteins(peptides, annot)
no_pairs <- data.table::data.table(V1 = character(), V2 = character())

test_that("with no alignments the score is the hit count", {
  hits <- data.table::data.table(u_pep_id = peptides,
                                 s1 = c(1, 1, 1, 1), s2 = c(1, 0, 0, 0))
  out <- polyclonal_scores(hits, pro, no_pairs)
  expect_equal(out$pro_id, c("P", "Q"))
  expect_equal(out$s1, c(3, 1))
  expect_equal(out$s2, c(1, 0))
})

test_that("two peptides that align collapse to one epitope", {
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(1, 1, 0, 0))
  pairs <- data.table::data.table(V1 = peptides[1L], V2 = peptides[2L])
  out <- polyclonal_scores(hits, pro, pairs)
  expect_equal(out[pro_id == "P", s1], 1)
})

test_that("an unaligned third hit is counted separately", {
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(1, 1, 1, 0))
  pairs <- data.table::data.table(V1 = peptides[1L], V2 = peptides[2L])
  out <- polyclonal_scores(hits, pro, pairs)
  # {p1,p2} align -> 1 epitope; p3 is independent -> 1. Total 2.
  expect_equal(out[pro_id == "P", s1], 2)
})

test_that("alignments only reduce the score when both peptides are hits", {
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(1, 0, 1, 0))
  pairs <- data.table::data.table(V1 = peptides[1L], V2 = peptides[2L])
  out <- polyclonal_scores(hits, pro, pairs)
  expect_equal(out[pro_id == "P", s1], 2)
})

test_that("samples with no hits score zero", {
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(0, 0, 0, 0))
  pairs <- data.table::data.table(V1 = peptides[1L], V2 = peptides[2L])
  out <- polyclonal_scores(hits, pro, pairs)
  expect_equal(out$s1, c(0, 0))
})

test_that("alignment pairs crossing proteins are ignored", {
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(1, 0, 0, 1))
  # p1 is in protein P, p4 in protein Q; the pair must not merge them.
  pairs <- data.table::data.table(V1 = peptides[1L], V2 = peptides[4L])
  out <- polyclonal_scores(hits, pro, pairs)
  expect_equal(out$s1, c(1, 1))
})

test_that("scores are per sample and memoisation does not leak between them", {
  hits <- data.table::data.table(u_pep_id = peptides,
                                 s1 = c(1, 1, 0, 0),
                                 s2 = c(1, 1, 1, 0),
                                 s3 = c(1, 1, 0, 0))
  pairs <- data.table::data.table(V1 = peptides[1L], V2 = peptides[2L])
  out <- polyclonal_scores(hits, pro, pairs)
  expect_equal(as.numeric(out[pro_id == "P", .(s1, s2, s3)]), c(1, 2, 1))
})

test_that("unannotated peptides are excluded", {
  partial <- annot[1:3]
  p <- suppressWarnings(map_peptides_to_proteins(peptides, partial))
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(1, 1, 1, 1))
  out <- polyclonal_scores(hits, p, no_pairs)
  expect_equal(out$pro_id, "P")
  expect_equal(out$s1, 3)
})

test_that("a mis-sized protein vector is rejected", {
  hits <- data.table::data.table(u_pep_id = peptides, s1 = c(1, 1, 1, 1))
  expect_error(polyclonal_scores(hits, pro[1:2], no_pairs), "length 2 but hits has 4")
})
