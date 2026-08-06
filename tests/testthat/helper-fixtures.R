# Shared fixtures: a miniature two-library screen. Some fixtures deliberately
# leave peptides out of the annotation, since that is the case that exposes
# peptide -> protein mapping errors.

ANNOT_FIELDS <- c("pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa",
                  "taxon_genus", "taxon_species", "gene_symbol", "product")

make_peptides <- function(lib, n, start = 1L) {
  sprintf("%s_%s", lib, formatC(seq(start, length.out = n), width = 6, flag = "0"))
}

# A matrix of `n` peptides x `n_samples` samples.
make_matrix <- function(peptides, n_samples = 3L, values = NULL, id = "u_pep_id") {
  dt <- data.table::data.table(x = peptides)
  data.table::setnames(dt, "x", id)
  for (s in seq_len(n_samples)) {
    v <- if (is.null(values)) as.numeric(seq_along(peptides) * s) else values[, s]
    data.table::set(dt, j = paste0("sample", s), value = v)
  }
  dt[]
}

# Annotation covering `peptides`, optionally omitting some of them.
make_annotation <- function(peptides, pro_id, omit = character()) {
  keep <- !(peptides %in% omit)
  dt <- data.table::data.table(u_pep_id = peptides[keep], pro_id = pro_id[keep])
  for (f in ANNOT_FIELDS)
    data.table::set(dt, j = f, value = paste0(f, "_", seq_len(sum(keep))))
  # A protein-level annotation table is also needed to annotate collapsed output.
  dt[]
}

# Write a metadata tree of the shape read_annotation() expects.
write_metadata <- function(root, libname, annot, pairs = NULL) {
  d <- file.path(root, libname)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(annot, file.path(d, paste0(libname, "_Universal.csv")))
  if (!is.null(pairs))
    data.table::fwrite(pairs, file.path(d, paste0(libname, "_trimmedpairs.csv")))
  d
}
