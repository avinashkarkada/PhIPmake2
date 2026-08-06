# Output comparison, for checking a run against a previous one.
#
# Differences are reported per file rather than treated as failures, since a
# protein-level table can legitimately change when the peptide -> protein
# mapping does.

#' Compare two phipmake output directories
#'
#' @param old Directory of reference output.
#' @param new Directory of output to check.
#' @param pattern Restrict to files matching this regex.
#' @param tolerance Absolute tolerance for numeric comparison.
#' @param max_report Maximum differing cells to record per file.
#' @return A `data.table`, one row per compared file, with the verdict and
#'   difference counts. Returned invisibly when `verbose = TRUE`.
#' @param verbose Print a summary.
#' @export
compare_outputs <- function(old, new, pattern = "\\.tsv$",
                            tolerance = 1e-8, max_report = 5L, verbose = TRUE) {
  stopifnot(dir.exists(old), dir.exists(new))
  rel <- function(d) {
    f <- list.files(d, pattern = pattern, recursive = TRUE)
    f[order(f)]
  }
  a <- rel(old); b <- rel(new)
  common <- intersect(a, b)

  rows <- lapply(common, function(f) {
    res <- .compare_one(file.path(old, f), file.path(new, f), tolerance, max_report)
    data.table::data.table(
      file = f, verdict = res$verdict, bytes = res$bytes, rows_old = res$rows_old,
      rows_new = res$rows_new, cols_old = res$cols_old, cols_new = res$cols_new,
      n_diff_cells = res$n_diff, max_abs_diff = res$max_abs, detail = res$detail)
  })

  out <- if (length(rows)) data.table::rbindlist(rows) else
    data.table::data.table(file = character(), verdict = character(),
                           bytes = character(),
                           rows_old = integer(), rows_new = integer(),
                           cols_old = integer(), cols_new = integer(),
                           n_diff_cells = integer(), max_abs_diff = numeric(),
                           detail = character())

  only_old <- setdiff(a, b); only_new <- setdiff(b, a)
  if (length(only_old))
    out <- rbind(out, data.table::data.table(
      file = only_old, verdict = "missing_in_new", bytes = NA_character_,
      rows_old = NA_integer_,
      rows_new = NA_integer_, cols_old = NA_integer_, cols_new = NA_integer_,
      n_diff_cells = NA_integer_, max_abs_diff = NA_real_, detail = ""))
  if (length(only_new))
    out <- rbind(out, data.table::data.table(
      file = only_new, verdict = "extra_in_new", bytes = NA_character_,
      rows_old = NA_integer_,
      rows_new = NA_integer_, cols_old = NA_integer_, cols_new = NA_integer_,
      n_diff_cells = NA_integer_, max_abs_diff = NA_real_, detail = ""))

  if (verbose) {
    tally <- table(out$verdict)
    message("Compared ", length(common), " file(s) present in both directories.")
    for (k in names(tally)) message(sprintf("  %-18s %d", k, tally[[k]]))
    nb <- sum(out$bytes == "same", na.rm = TRUE)
    message(sprintf("  %-18s %d of %d", "byte-identical", nb, length(common)))
    same_vals <- out[verdict == "identical" & bytes == "differ"]
    if (nrow(same_vals))
      message(sprintf("  %-18s %d (same values, different formatting)",
                      "formatting only", nrow(same_vals)))
    bad <- out[verdict %chin% c("values_differ", "shape_differs", "ids_differ")]
    if (nrow(bad))
      for (i in seq_len(min(nrow(bad), 10L)))
        message(sprintf("  ! %s: %s", bad$file[i], bad$detail[i]))
    return(invisible(out))
  }
  out
}

.compare_one <- function(pa, pb, tolerance, max_report) {
  blank <- list(verdict = NA_character_, bytes = NA_character_,
                rows_old = NA_integer_, rows_new = NA_integer_,
                cols_old = NA_integer_, cols_new = NA_integer_, n_diff = NA_integer_,
                max_abs = NA_real_, detail = "")

  # Byte comparison first. Parsing strips quoting and normalises number
  # formatting, so two files can hold identical values and still differ on disk.
  blank$bytes <- if (identical(file.size(pa), file.size(pb)) &&
                     identical(unname(tools::md5sum(pa)), unname(tools::md5sum(pb))))
    "same" else "differ"
  A <- try(data.table::fread(pa, showProgress = FALSE), silent = TRUE)
  B <- try(data.table::fread(pb, showProgress = FALSE), silent = TRUE)
  if (inherits(A, "try-error") || inherits(B, "try-error")) {
    blank$verdict <- "unreadable"; return(blank)
  }

  blank$rows_old <- nrow(A); blank$rows_new <- nrow(B)
  blank$cols_old <- ncol(A); blank$cols_new <- ncol(B)

  if (ncol(A) != ncol(B) || !identical(names(A), names(B))) {
    blank$verdict <- "shape_differs"
    blank$detail <- sprintf("columns %d vs %d", ncol(A), ncol(B))
    return(blank)
  }
  if (nrow(A) != nrow(B)) {
    blank$verdict <- "shape_differs"
    blank$detail <- sprintf("rows %d vs %d", nrow(A), nrow(B))
    return(blank)
  }
  if (!identical(as.character(A[[1L]]), as.character(B[[1L]]))) {
    n <- sum(as.character(A[[1L]]) != as.character(B[[1L]]))
    blank$verdict <- "ids_differ"
    blank$detail <- sprintf("%d identifier(s) differ or are reordered", n)
    return(blank)
  }

  n_diff <- 0L; max_abs <- 0; examples <- character()
  for (j in names(A)[-1L]) {
    x <- A[[j]]; y <- B[[j]]
    if (is.numeric(x) && is.numeric(y)) {
      d <- abs(x - y)
      bad <- which(!is.na(d) & d > tolerance | xor(is.na(x), is.na(y)))
    } else {
      bad <- which(as.character(x) != as.character(y) | xor(is.na(x), is.na(y)))
      d <- NULL
    }
    if (length(bad)) {
      n_diff <- n_diff + length(bad)
      if (!is.null(d)) max_abs <- max(max_abs, max(d[bad], na.rm = TRUE))
      if (length(examples) < max_report)
        examples <- c(examples, sprintf("%s[%s]: %s vs %s", j, A[[1L]][bad[1L]],
                                        format(x[bad[1L]]), format(y[bad[1L]])))
    }
  }

  blank$n_diff <- n_diff
  blank$max_abs <- max_abs
  if (n_diff == 0L) {
    blank$verdict <- "identical"
  } else {
    blank$verdict <- "values_differ"
    blank$detail <- sprintf("%d cell(s); %s", n_diff,
                            paste(utils::head(examples, max_report), collapse = "; "))
  }
  blank
}
