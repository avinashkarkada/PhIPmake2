# Joining peptide / protein annotations onto a data matrix.

.PEPTIDE_FIELDS <- c("pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa",
                     "taxon_genus", "taxon_species", "gene_symbol", "product")
# Fields that only make sense per peptide; dropped for protein-level tables.
.PEPTIDE_ONLY_FIELDS <- c("pep_id", "pos_start", "pos_end", "pep_aa")

#' Attach annotation columns to a data matrix
#'
#' Column order is identifier, annotation fields, then samples. Annotation
#' values are written as text; because this runs one sub-library at a time, the
#' character columns never exist for more than one library at once.
#'
#' @param data Matrix whose first column is `u_pep_id` or `pro_id`.
#' @param annot Annotation table for the same library.
#' @param fields Annotation columns to attach.
#' @return A `data.table`.
#' @export
annotate_table <- function(data, annot, fields = .PEPTIDE_FIELDS) {
  data <- data.table::as.data.table(data)
  annot <- data.table::as.data.table(annot)

  id_col <- names(data)[1L]
  if (identical(id_col, "pro_id")) fields <- setdiff(fields, .PEPTIDE_ONLY_FIELDS)

  if (!id_col %in% names(annot))
    stop("Annotation has no '", id_col, "' column to join on; found: ",
         paste(utils::head(names(annot), 12L), collapse = ", "), call. = FALSE)

  missing_fields <- setdiff(fields, names(annot))
  if (length(missing_fields))
    stop("Annotation is missing requested field(s): ",
         paste(missing_fields, collapse = ", "), call. = FALSE)

  i <- match(as.character(data[[1L]]), as.character(annot[[id_col]]))
  found <- !is.na(i)

  out <- data.table::data.table(x = data[[1L]])
  data.table::setnames(out, "x", id_col)
  for (f in fields) {
    v <- as.character(annot[[f]][i])
    # A value that is missing in the annotation is written as the literal "NA",
    # which keeps the column readable as text. A peptide with no annotation row
    # at all stays a true NA, so the two cases remain distinguishable.
    v[found & is.na(v)] <- "NA"
    data.table::set(out, j = f, value = v)
  }
  for (v in names(data)[-1L])
    data.table::set(out, j = v, value = data[[v]])
  out[]
}
