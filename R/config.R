# Run configuration.
#
# A plain list. The drake_params.tsv that screen directories already carry is
# read as-is, so existing directories work without modification.

.DEFAULT_PARAMS <- list(
  counts_filename     = "counts.csv",
  counts_type         = "Counts",
  foldchange_filename = "fold_change.csv",
  foldchange_type     = "FoldChange",
  enrichment_filename = "enrichment.csv",
  enrichment_type     = "EdgeR",
  hits_filename       = "Hits.csv",
  metadata_path       = Sys.getenv("PHIPMAKE2_METADATA", unset = NA_character_),
  output_extension    = "tsv"
)

#' Build a phipmake2 run configuration
#'
#' @param wd Screen directory holding the input matrices.
#' @param screen_name Output prefix. Defaults to `basename(wd)`.
#' @param params_path Optional `drake_params.tsv` to read defaults from. When
#'   NULL, `wd/drake_params.tsv` is used if it exists.
#' @param ... Any parameter in `.DEFAULT_PARAMS`, overriding the file.
#' @return A list of resolved settings.
#' @export
phipmake2_config <- function(wd, screen_name = NULL, params_path = NULL, ...) {
  wd <- normalizePath(wd, mustWork = TRUE)
  cfg <- .DEFAULT_PARAMS

  if (is.null(params_path)) {
    candidate <- file.path(wd, "drake_params.tsv")
    if (file.exists(candidate)) params_path <- candidate
  }
  if (!is.null(params_path) && file.exists(params_path)) {
    p <- data.table::fread(params_path, showProgress = FALSE)
    kv <- stats::setNames(as.character(p[[2L]]), as.character(p[[1L]]))
    for (k in intersect(names(kv), c(names(cfg), "screen_name")))
      cfg[[k]] <- kv[[k]]
    if (is.null(screen_name) && !is.na(kv["screen_name"]))
      screen_name <- unname(kv["screen_name"])
  }

  overrides <- list(...)
  unknown <- setdiff(names(overrides), c(names(.DEFAULT_PARAMS), "screen_name"))
  if (length(unknown))
    stop("Unknown configuration parameter(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  for (k in names(overrides)) cfg[[k]] <- overrides[[k]]

  cfg$wd <- wd
  cfg$screen_name <- if (!is.null(screen_name)) screen_name else basename(wd)
  cfg
}

# Absolute path of an input matrix, or NA when absent.
.input_path <- function(cfg, which) {
  p <- file.path(cfg$wd, cfg[[paste0(which, "_filename")]])
  if (file.exists(p)) p else NA_character_
}

# Output path builder. `stem` is the part between the screen name and the
# extension, e.g. "Counts", "Hits_foldchange_promax".
.out_path <- function(cfg, stem, lib = NULL, annotated = FALSE) {
  suffix <- if (annotated) "_annotated" else ""
  ext <- paste0(".", cfg$output_extension)
  if (is.null(lib)) {
    file.path(cfg$wd, paste0(cfg$screen_name, "_", stem, suffix, ext))
  } else {
    file.path(cfg$wd, lib,
              paste0(cfg$screen_name, "_", lib, "_", stem, suffix, ext))
  }
}
