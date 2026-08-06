# The pipeline.
#
# Peak memory is what drives the design here. A screen has four input matrices
# and produces a dozen derived tables from each; holding those simultaneously is
# what makes large screens fail. So:
#
#   * one input matrix is loaded at a time and freed before the next;
#   * within a matrix, one sub-library is processed at a time;
#   * whole-screen ("pan") files are appended to disk per sub-library rather
#     than assembled in memory;
#   * nothing is cached. Re-running is cheap because the work is cheap.
#
# Hits are the exception: the masked variants of counts, fold change and
# enrichment all need them, so the hit matrix stays loaded across those passes.

.info <- function(verbose, fmt, ...) if (verbose) message(sprintf(fmt, ...))

# Recognised in a plan string but not implemented. AVARDA is the one we expect
# to add; epitopefindr and pairwise were never functional upstream.
.PLANNED_STAGES <- c("avarda", "epitopefindr", "pairwise")

# A pair of pan writers (raw + annotated) for one output family.
.stream <- function(cfg, stem) {
  list(stem = stem,
       raw = .pan_writer(.out_path(cfg, stem)),
       ann = .pan_writer(.out_path(cfg, stem, annotated = TRUE)))
}

# Write one sub-library table: per-library file plus append to both pan files.
.stream_add <- function(st, cfg, tbl, lib, annot, pan_raw = TRUE) {
  write_matrix(tbl, .out_path(cfg, st$stem, lib = lib))
  if (pan_raw) st$raw$add(tbl)
  ann <- annotate_table(tbl, annot)
  write_matrix(ann, .out_path(cfg, st$stem, lib = lib, annotated = TRUE))
  st$ann$add(ann)
  invisible(NULL)
}

# Load an input matrix, write its verbatim pan copy, and normalise the id column.
.load_input <- function(cfg, which, stem, verbose) {
  path <- .input_path(cfg, which)
  if (is.na(path)) return(NULL)
  .info(verbose, "Reading %s", basename(path))
  dt <- read_matrix(path, id_col = NULL)
  # The pan file for a base data type is the input re-serialised under its
  # original first-column name, not a concatenation of sub-libraries.
  write_matrix(dt, .out_path(cfg, stem))
  data.table::setnames(dt, 1L, "u_pep_id")
  dt
}

#' Run the phipmake2 pipeline
#'
#' @param wd Screen directory containing the input matrices.
#' @param stages Which data types to process. Any of `"counts"`,
#'   `"foldchange"`, `"enrichment"`, `"hits"`, `"polyclonal"`. Also accepts a
#'   single dash-separated string as used by the original pipeline, e.g.
#'   `"Counts-FoldChange-Enrichment-Hits"`.
#' @param screen_name Output prefix; defaults to `basename(wd)`.
#' @param metadata_path Directory of per-library annotation folders.
#' @param threads `data.table` threads; defaults to the current setting.
#' @param verbose Print progress.
#' @param ... Further overrides passed to [phipmake2_config()].
#' @return Invisibly, a character vector of files written.
#' @export
run_phipmake2 <- function(wd,
                          stages = c("counts", "foldchange", "enrichment", "hits"),
                          screen_name = NULL,
                          metadata_path = NULL,
                          threads = NULL,
                          verbose = TRUE,
                          ...) {
  if (length(stages) == 1L && grepl("-", stages, fixed = TRUE))
    stages <- strsplit(stages, "-", fixed = TRUE)[[1L]]
  stages <- tolower(stages)
  valid <- c("counts", "foldchange", "enrichment", "hits", "polyclonal")

  # Stages the original pipeline accepted that are not implemented yet. They are
  # skipped with a warning rather than rejected, so existing submission scripts
  # carrying them keep working.
  planned <- intersect(stages, .PLANNED_STAGES)
  if (length(planned)) {
    warning("Stage(s) not implemented and skipped: ", paste(planned, collapse = ", "),
            ".", call. = FALSE)
    stages <- setdiff(stages, planned)
  }

  bad <- setdiff(stages, valid)
  if (length(bad))
    stop("Unknown stage(s): ", paste(bad, collapse = ", "),
         ". Valid stages: ", paste(valid, collapse = ", "), call. = FALSE)
  if (length(stages) == 0L)
    stop("No runnable stages requested.", call. = FALSE)

  if (!is.null(threads)) data.table::setDTthreads(threads)

  extra <- list(...)
  if (!is.null(metadata_path)) extra$metadata_path <- metadata_path
  cfg <- do.call(phipmake2_config,
                 c(list(wd = wd, screen_name = screen_name), extra))

  t0 <- Sys.time()
  .info(verbose, "phipmake2 %s | screen '%s'",
        utils::packageVersion("phipmake2"), cfg$screen_name)
  .info(verbose, "Stages: %s", paste(stages, collapse = ", "))

  # ---- Establish sub-libraries from the identifier column only ----------
  ref <- NULL
  for (w in c("counts", "enrichment", "foldchange", "hits")) {
    p <- .input_path(cfg, w)
    if (!is.na(p)) { ref <- p; break }
  }
  if (is.null(ref))
    stop("No input matrices found in ", cfg$wd,
         ". Expected one of: ", paste(vapply(c("counts", "foldchange", "enrichment", "hits"),
                                             function(w) cfg[[paste0(w, "_filename")]], ""),
                                      collapse = ", "), call. = FALSE)

  ids <- data.table::fread(ref, select = 1L, showProgress = FALSE)[[1L]]
  lib_index <- peptide_library_index(ids)
  libs <- names(lib_index)
  .info(verbose, "Libraries (%d): %s", length(libs), paste(libs, collapse = ", "))
  for (l in libs) dir.create(file.path(cfg$wd, l), showWarnings = FALSE, recursive = TRUE)

  .info(verbose, "Loading annotations from %s", cfg$metadata_path)
  annots <- lapply(libs, read_annotation, metadata_path = cfg$metadata_path)
  names(annots) <- libs
  pro_maps <- lapply(libs, function(l)
    map_peptides_to_proteins(ids[lib_index[[l]]], annots[[l]]))
  names(pro_maps) <- libs

  slice <- function(dt, l) dt[lib_index[[l]]]

  # ---- Hits are needed by the masked variants of every other type -------
  hits <- NULL
  want_hits <- "hits" %in% stages || "polyclonal" %in% stages
  if (want_hits) {
    hits <- .load_input(cfg, "hits", "Hits", verbose)
    if (is.null(hits))
      stop("Stage 'hits' requested but ", cfg$hits_filename, " is not in ", cfg$wd,
           call. = FALSE)
  }

  # ---- counts -----------------------------------------------------------
  if ("counts" %in% stages) {
    counts <- .load_input(cfg, "counts", cfg$counts_type, verbose)
    if (is.null(counts)) {
      warning("Stage 'counts' skipped: ", cfg$counts_filename, " not found.", call. = FALSE)
    } else {
      s_raw    <- .stream(cfg, cfg$counts_type)
      s_prosum <- .stream(cfg, paste0(cfg$counts_type, "_prosum"))
      s_hitc   <- if (!is.null(hits)) .stream(cfg, "Hits_counts") else NULL
      for (l in libs) {
        .info(verbose, "counts | %s", l)
        sub <- slice(counts, l)
        .stream_add(s_raw, cfg, sub, l, annots[[l]], pan_raw = FALSE)
        .stream_add(s_prosum, cfg,
                    collapse_peptides(sub, pro_maps[[l]], "sum"), l, annots[[l]])
        if (!is.null(s_hitc))
          .stream_add(s_hitc, cfg, emphasize_hits(sub, slice(hits, l), 0), l, annots[[l]])
        rm(sub); invisible(gc(FALSE))
      }
      rm(counts); invisible(gc(FALSE))
    }
  }

  # ---- fold change ------------------------------------------------------
  if ("foldchange" %in% stages) {
    fc <- .load_input(cfg, "foldchange", cfg$foldchange_type, verbose)
    if (is.null(fc)) {
      warning("Stage 'foldchange' skipped: ", cfg$foldchange_filename, " not found.", call. = FALSE)
    } else {
      s_raw  <- .stream(cfg, cfg$foldchange_type)
      s_hfc  <- if (!is.null(hits)) .stream(cfg, "Hits_foldchange") else NULL
      s_hfcm <- if (!is.null(hits)) .stream(cfg, "Hits_foldchange_promax") else NULL
      for (l in libs) {
        .info(verbose, "foldchange | %s", l)
        sub <- slice(fc, l)
        .stream_add(s_raw, cfg, sub, l, annots[[l]], pan_raw = FALSE)
        if (!is.null(s_hfc)) {
          masked <- emphasize_hits(sub, slice(hits, l), 1)
          .stream_add(s_hfc, cfg, masked, l, annots[[l]])
          .stream_add(s_hfcm, cfg,
                      collapse_peptides(masked, pro_maps[[l]], "max"), l, annots[[l]])
          rm(masked)
        }
        rm(sub); invisible(gc(FALSE))
      }
      rm(fc); invisible(gc(FALSE))
    }
  }

  # ---- enrichment -------------------------------------------------------
  if ("enrichment" %in% stages) {
    enr <- .load_input(cfg, "enrichment", cfg$enrichment_type, verbose)
    if (is.null(enr)) {
      warning("Stage 'enrichment' skipped: ", cfg$enrichment_filename, " not found.", call. = FALSE)
    } else {
      s_raw  <- .stream(cfg, cfg$enrichment_type)
      s_pmax <- .stream(cfg, paste0(cfg$enrichment_type, "_promax"))
      s_he   <- if (!is.null(hits)) .stream(cfg, "Hits_enrichment") else NULL
      s_hem  <- if (!is.null(hits)) .stream(cfg, "Hits_enrichment_promax") else NULL
      for (l in libs) {
        .info(verbose, "enrichment | %s", l)
        sub <- slice(enr, l)
        .stream_add(s_raw, cfg, sub, l, annots[[l]], pan_raw = FALSE)
        .stream_add(s_pmax, cfg,
                    collapse_peptides(sub, pro_maps[[l]], "max"), l, annots[[l]])
        if (!is.null(s_he)) {
          masked <- emphasize_hits(sub, slice(hits, l), 0)
          .stream_add(s_he, cfg, masked, l, annots[[l]])
          .stream_add(s_hem, cfg,
                      collapse_peptides(masked, pro_maps[[l]], "max"), l, annots[[l]])
          rm(masked)
        }
        rm(sub); invisible(gc(FALSE))
      }
      rm(enr); invisible(gc(FALSE))
    }
  }

  # ---- hits and polyclonal ---------------------------------------------
  if ("hits" %in% stages) {
    s_raw    <- .stream(cfg, "Hits")
    s_prosum <- .stream(cfg, "Hits_prosum")
    for (l in libs) {
      .info(verbose, "hits | %s", l)
      sub <- slice(hits, l)
      .stream_add(s_raw, cfg, sub, l, annots[[l]], pan_raw = FALSE)
      .stream_add(s_prosum, cfg,
                  collapse_peptides(sub, pro_maps[[l]], "sum"), l, annots[[l]])
      rm(sub); invisible(gc(FALSE))
    }
  }

  if ("polyclonal" %in% stages) {
    s_poly <- .stream(cfg, "Polyclonal")
    for (l in libs) {
      .info(verbose, "polyclonal | %s", l)
      pairs <- read_pairs(l, cfg$metadata_path)
      if (nrow(pairs) == 0L)
        warning("No alignment pairs for '", l,
                "'; polyclonal score reduces to the hit count.", call. = FALSE)
      sub <- slice(hits, l)
      .stream_add(s_poly, cfg,
                  polyclonal_scores(sub, pro_maps[[l]], pairs, verbose = verbose),
                  l, annots[[l]])
      rm(sub, pairs); invisible(gc(FALSE))
    }
  }

  elapsed <- difftime(Sys.time(), t0, units = "mins")
  .info(verbose, "Done in %.1f min", as.numeric(elapsed))
  invisible(list.files(cfg$wd, pattern = paste0("^", cfg$screen_name),
                       recursive = TRUE, full.names = TRUE))
}
