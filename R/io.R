# Reading and writing of PhIP-Seq matrices and metadata.
#
# Whole-screen ("pan") outputs are never held in memory. Each sub-library is
# appended to the pan file as it is produced, so the pan file costs nothing
# beyond the sub-library already in hand.

#' Read a PhIP-Seq matrix
#'
#' @param path Path to a counts / fold change / enrichment / hits matrix.
#' @param id_col Name to coerce the first (identifier) column to, or NULL to
#'   preserve whatever the file declares.
#' @param nThread Threads for `data.table::fread`.
#' @return A `data.table`.
#' @export
read_matrix <- function(path, id_col = "u_pep_id", nThread = data.table::getDTthreads()) {
  if (!file.exists(path)) stop("Matrix file not found: ", path, call. = FALSE)
  dt <- data.table::fread(path, showProgress = FALSE, nThread = nThread)
  if (!is.null(id_col)) data.table::setnames(dt, 1L, id_col)
  dt
}

#' Write a PhIP-Seq matrix
#'
#' Tab separated, literal `NA` for missing, column names and character values
#' quoted, numeric values bare.
#'
#' `quote = TRUE` rather than the "auto" default: recent `data.table` leaves
#' column names unquoted under "auto" where older versions quoted them, and the
#' header is part of the file format these outputs are expected to have.
#'
#' @param x A `data.table` or `data.frame`.
#' @param path Destination path. Parent directories are created as needed.
#' @param append If TRUE, append without repeating the header.
#' @return `path`, invisibly.
#' @export
write_matrix <- function(x, path, append = FALSE) {
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "NA", quote = TRUE,
                     append = append, col.names = !append,
                     showProgress = FALSE)
  invisible(path)
}

# Internal: a small handle that turns a sequence of per-library writes into one
# pan file, without ever concatenating in memory.
.pan_writer <- function(path) {
  started <- FALSE
  list(
    add = function(x) {
      write_matrix(x, path, append = started)
      started <<- TRUE
      invisible(NULL)
    },
    started = function() started
  )
}

#' Read a peptide annotation file for one library
#'
#' The library directory is matched exactly; a missing directory is an error
#' rather than a fallback to whatever else happens to be there.
#'
#' @param libname Resolved library name, e.g. `LibAlpha_000`.
#' @param metadata_path Directory holding one sub-directory per library.
#' @param na_as_string If TRUE (default) missing annotation values become the
#'   literal string "NA".
#' @return A `data.table` of annotations.
#' @export
read_annotation <- function(libname, metadata_path, na_as_string = TRUE) {
  dir <- .resolve_library_dir(libname, metadata_path)
  hit <- .pick_library_file(dir, libname, "_Universal")
  if (length(hit) == 0L)
    stop("No '_Universal' annotation file for library '", libname, "' in ", dir, call. = FALSE)

  annot <- data.table::fread(hit, showProgress = FALSE)
  # Annotation values are emitted as text downstream, so a missing value is
  # carried through as the literal string rather than a typed NA.
  if (na_as_string) {
    for (j in names(annot)) {
      col <- annot[[j]]
      if (anyNA(col)) data.table::set(annot, which(is.na(col)), j,
                                      if (is.character(col)) "NA" else NA)
    }
  }
  annot
}

#' Read intra-library BLAST alignment pairs used for polyclonal scoring
#'
#' @param libname Resolved library name.
#' @param pairs_path Directory holding one sub-directory per library.
#' @param pattern Substring identifying the pairs file.
#' @return A two-column `data.table` of aligning peptide pairs, or an empty
#'   table when the library has no pairs file.
#' @export
read_pairs <- function(libname, pairs_path, pattern = "_trimmedpairs") {
  dir <- tryCatch(.resolve_library_dir(libname, pairs_path), error = function(e) NA_character_)
  empty <- data.table::data.table(V1 = character(), V2 = character())
  if (is.na(dir)) return(empty)
  hit <- .pick_library_file(dir, libname, pattern)
  if (length(hit) == 0L) return(empty)
  data.table::fread(hit, showProgress = FALSE)
}

.base_of <- function(libname) sub("_[0-9]+$", "", libname)

# Pick one file out of a library's metadata directory.
#
# Filenames are not consistent across libraries: most carry a date prefix
# (20190530_LibAlpha_001_Universal.tsv) and some do not
# (LibTheta_001_Universal.tsv), so the library name cannot be assumed to
# be at the front. The directory has already been resolved exactly, so `kind`
# ("_Universal", "_trimmedpairs") is usually enough on its own. Where a
# directory holds more than one candidate, one naming this exact library version
# wins; failing that the last in sort order does, which takes the most recent
# date prefix.
.pick_library_file <- function(dir, libname, kind) {
  files <- list.files(dir, full.names = TRUE)
  files <- files[!dir.exists(files)]           # skip Historic/ and friends
  hit <- files[grepl(kind, basename(files), fixed = TRUE)]
  if (length(hit) <= 1L) return(hit)

  exact <- hit[grepl(libname, basename(hit), fixed = TRUE)]
  if (length(exact) == 1L) return(exact)
  if (length(exact) > 1L) hit <- exact

  chosen <- hit[order(basename(hit))][length(hit)]
  warning("Multiple '", kind, "' files for library '", libname, "' in ", dir,
          "; using ", basename(chosen), ".", call. = FALSE)
  chosen
}

# Locate a library's metadata directory by exact name, falling back to the
# base name (a merged LibAlpha_000 is documented under LibAlpha_001/2).
.resolve_library_dir <- function(libname, root) {
  if (!dir.exists(root)) stop("Metadata directory not found: ", root, call. = FALSE)
  dirs <- list.dirs(root, full.names = FALSE, recursive = FALSE)
  if (libname %chin% dirs) return(file.path(root, libname))

  base <- .base_of(libname)
  cand <- dirs[startsWith(dirs, paste0(base, "_")) | dirs == base]
  if (length(cand) == 0L)
    stop("No metadata directory for library '", libname, "' under ", root,
         ". Available: ", paste(utils::head(dirs, 10L), collapse = ", "), call. = FALSE)
  file.path(root, cand[order(cand)][1L])
}
