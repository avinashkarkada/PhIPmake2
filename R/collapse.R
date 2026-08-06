# Peptide -> protein collapsing.
#
# prosum and promax both reduce a peptide x sample matrix to a protein x sample
# matrix and differ only in the aggregation function, so they share one grouped
# data.table pass.
#
# The subtlety worth knowing about is the peptide -> protein map: it has to stay
# the same length as the data, carrying NA for peptides with no annotation.
# Dropping those NAs first gives a shorter vector, and R will happily recycle it
# when it is used as a grouping vector -- which quietly assigns peptides to the
# wrong proteins instead of raising an error.

#' Map peptides to proteins
#'
#' Returns a vector the same length as `u_pep_id`, with `NA` for peptides that
#' are absent from the annotation, so it can be used as a grouping vector
#' without risk of recycling.
#'
#' @param u_pep_id Character vector of peptide identifiers, in data row order.
#' @param annot Annotation table containing `u_pep_id` and `pro_id`.
#' @param warn_unmapped Warn with a count of unmapped peptides.
#' @return Character vector of protein identifiers, parallel to `u_pep_id`.
#' @export
map_peptides_to_proteins <- function(u_pep_id, annot, warn_unmapped = TRUE) {
  if (!all(c("u_pep_id", "pro_id") %in% names(annot)))
    stop("Annotation must contain 'u_pep_id' and 'pro_id' columns; found: ",
         paste(names(annot), collapse = ", "), call. = FALSE)

  pro <- as.character(annot$pro_id)[match(as.character(u_pep_id),
                                          as.character(annot$u_pep_id))]
  n_missing <- sum(is.na(pro))
  if (warn_unmapped && n_missing > 0L)
    warning(sprintf(
      "%d of %d peptides have no annotation and are excluded from protein-level output.",
      n_missing, length(pro)), call. = FALSE)
  pro
}

#' Collapse a peptide matrix to a protein matrix
#'
#' @param data Peptide matrix; first column identifiers, remaining columns samples.
#' @param pro_id Protein identifier per row, as returned by
#'   [map_peptides_to_proteins()]. Must be `nrow(data)` long.
#' @param fun `"sum"` for prosum or `"max"` for promax.
#' @param na_rows `"drop"` (default) discards any peptide with a missing value in
#'   any sample before aggregating. `"keep"` aggregates per sample with
#'   `na.rm = TRUE` instead.
#' @return A `data.table` with a `pro_id` column followed by the sample columns.
#' @export
collapse_peptides <- function(data, pro_id, fun = c("sum", "max"),
                              na_rows = c("drop", "keep")) {
  fun <- match.arg(fun)
  na_rows <- match.arg(na_rows)
  data <- data.table::as.data.table(data)

  if (length(pro_id) != nrow(data))
    stop(sprintf(
      "pro_id has length %d but data has %d rows. Pass the full-length vector from map_peptides_to_proteins() rather than dropping NAs from it first.",
      length(pro_id), nrow(data)), call. = FALSE)

  value_cols <- names(data)[-1L]
  if (length(value_cols) == 0L)
    stop("data has no sample columns to collapse.", call. = FALSE)

  mapped <- !is.na(pro_id)
  # Proteins appear in order of first appearance among annotated peptides.
  proteins <- unique(pro_id[mapped])
  if (length(proteins) == 0L) {
    out <- data.table::data.table(pro_id = character())
    for (j in value_cols) data.table::set(out, j = j, value = numeric())
    return(out[])
  }

  # Which rows contribute values. Tested column by column so we never
  # materialise a copy of the sample block just to look for NAs.
  contributes <- mapped
  if (na_rows == "drop") {
    for (j in value_cols) {
      col <- data[[j]]
      if (anyNA(col)) contributes <- contributes & !is.na(col)
    }
  }

  grp <- ifelse(contributes, pro_id, NA_character_)
  agg <- if (fun == "sum") function(x) sum(x, na.rm = TRUE) else function(x) max(x, na.rm = TRUE)

  # Grouping on an external vector keeps `data` from being copied.
  out <- data[, lapply(.SD, agg), by = list(pro_id = grp), .SDcols = value_cols]

  # Re-index onto the full protein set so that a protein whose every peptide was
  # dropped comes back as NA rather than disappearing from the output.
  out <- out[match(proteins, out$pro_id)]
  data.table::set(out, j = "pro_id", value = proteins)
  out[]
}

#' Replace non-hit values with a constant
#'
#' @param data Matrix of values; first column identifiers.
#' @param hits Binary hit matrix aligned to `data`, same shape and row order.
#' @param default Value written wherever `hits` is 0.
#' @return A `data.table` of the same shape as `data`.
#' @export
emphasize_hits <- function(data, hits, default) {
  data <- data.table::as.data.table(data)
  hits <- data.table::as.data.table(hits)

  if (nrow(data) != nrow(hits))
    stop(sprintf("data has %d rows but hits has %d.", nrow(data), nrow(hits)), call. = FALSE)
  if (ncol(data) != ncol(hits))
    stop(sprintf("data has %d columns but hits has %d.", ncol(data), ncol(hits)), call. = FALSE)
  if (!identical(as.character(data[[1L]]), as.character(hits[[1L]])))
    stop("data and hits are not in the same row order; refusing to mask by position.",
         call. = FALSE)

  out <- data.table::copy(data)
  value_cols <- names(out)[-1L]
  hit_cols <- names(hits)[-1L]
  for (k in seq_along(value_cols)) {
    zero <- which(hits[[hit_cols[k]]] == 0)
    if (length(zero)) data.table::set(out, i = zero, j = value_cols[k], value = default)
  }
  out[]
}
