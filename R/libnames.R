# Parsing of universal peptide identifiers (u_pep_id).
#
# A u_pep_id looks like "LibAlpha_001_ATGCGC...": a library base name, an
# underscore, a zero-padded numeric library version, an underscore, and the
# nucleotide sequence.
#
# Splitting on "_" and assuming three fields does not hold in general, since
# both the base name and the sequence can themselves contain underscores.

# Anchored, lazy on the base so that base names containing underscores still
# parse: the version is the first all-digit field, everything before it is the
# base, everything after is the sequence.
.U_PEP_ID_RE <- "^(.+?)_([0-9]+)_(.+)$"

#' Split universal peptide identifiers into base, version and sequence
#'
#' @param u_pep_id Character vector of universal peptide identifiers.
#' @param strict If TRUE (default) unparseable identifiers raise an error rather
#'   than being returned as NA.
#' @return A `data.table` with columns `u_pep_id`, `base`, `version`, `sequence`
#'   and `library` (the versioned name, e.g. `LibAlpha_001`).
#' @export
parse_u_pep_id <- function(u_pep_id, strict = TRUE) {
  u_pep_id <- as.character(u_pep_id)
  m <- regmatches(u_pep_id, regexec(.U_PEP_ID_RE, u_pep_id))
  ok <- lengths(m) == 4L

  if (!all(ok)) {
    bad <- utils::head(u_pep_id[!ok], 5L)
    msg <- sprintf(
      "%d of %d peptide identifiers do not match '<base>_<version>_<sequence>'. First offenders: %s",
      sum(!ok), length(ok), paste(bad, collapse = ", "))
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  base <- version <- sequence <- rep(NA_character_, length(u_pep_id))
  if (any(ok)) {
    got <- do.call(rbind, m[ok])
    base[ok]     <- got[, 2L]
    version[ok]  <- got[, 3L]
    sequence[ok] <- got[, 4L]
  }

  data.table::data.table(
    u_pep_id = u_pep_id,
    base     = base,
    version  = version,
    sequence = sequence,
    library  = ifelse(is.na(base), NA_character_, paste(base, version, sep = "_"))
  )
}

#' Resolve the output library name for each peptide
#'
#' A base library present at a single version keeps that version
#' (`LibGamma_001`). A base library present at two or more versions is merged
#' under version `000` (`LibAlpha_001` + `LibAlpha_002` -> `LibAlpha_000`).
#'
#' @param u_pep_id Character vector of universal peptide identifiers.
#' @return Character vector of resolved library names, parallel to `u_pep_id`.
#' @export
library_of <- function(u_pep_id) {
  p <- parse_u_pep_id(u_pep_id)
  n_versions <- p[, list(n = data.table::uniqueN(version)), by = base]
  merged <- n_versions[n > 1L, base]
  ifelse(p$base %chin% merged, paste0(p$base, "_000"), p$library)
}

#' Index peptides by resolved sub-library
#'
#' Parses identifiers once and returns row indices, so callers can subset each
#' data file without rescanning.
#'
#' @param u_pep_id Character vector of universal peptide identifiers.
#' @return A named list of integer row-index vectors, one per sub-library, in
#'   order of first appearance.
#' @export
peptide_library_index <- function(u_pep_id) {
  lib <- library_of(u_pep_id)
  idx <- split(seq_along(lib), factor(lib, levels = unique(lib)))
  idx[order(match(names(idx), unique(lib)))]
}
