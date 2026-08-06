# Polyclonal scoring: the minimum number of independent epitopes explaining a
# sample's hits within one protein.
#
# The scoring rule is Daniel Monaco's independence filter, unchanged. Computing
# it naively means one graph per (protein x sample) pair, which for a
# 20,000-protein library across 145 samples is ~2.9 million graph builds. Four
# things keep that manageable:
#
#   * proteins with no internal alignments skip the graph entirely -- their
#     score is just the hit count, and they are the large majority;
#   * each protein's alignment graph is built once and induced per sample;
#   * samples whose hits induce no edges short-circuit to the hit count;
#   * repeated hit sets within a protein are memoised.

.independence_number <- function(g) {
  if (igraph::vcount(g) == 0L) return(0L)
  as.integer(igraph::ivs_size(g))
}

# Faithful port of independence_filter(): components under 30 vertices are
# solved exactly; larger ones have their highest-degree vertex removed
# repeatedly until no vertex exceeds degree 5, then are solved per component.
.independent_hits <- function(g) {
  if (igraph::vcount(g) == 0L) return(0L)
  if (igraph::ecount(g) == 0L) return(igraph::vcount(g))

  if (igraph::is_directed(g)) g <- igraph::as_undirected(g, mode = "collapse")
  comps <- igraph::decompose(g)
  # vcount() returns a double in igraph 2.x, so do not demand an integer here.
  sizes <- vapply(comps, function(x) as.numeric(igraph::vcount(x)), numeric(1L))

  total <- 0L
  for (cmp in comps[sizes < 30L]) total <- total + .independence_number(cmp)

  for (cmp in comps[sizes >= 30L]) {
    deg <- igraph::degree(cmp)
    while (length(deg) && max(deg) > 5) {
      drop <- which.max(deg)
      cmp <- igraph::delete_vertices(cmp, igraph::V(cmp)[drop])
      deg <- igraph::degree(cmp)
    }
    for (sub in igraph::decompose(cmp)) total <- total + .independence_number(sub)
  }
  total
}

#' Polyclonal scores for one library
#'
#' @param hits Binary hit matrix; first column `u_pep_id`, remaining columns samples.
#' @param pro_id Protein per row, from [map_peptides_to_proteins()]; length `nrow(hits)`.
#' @param pairs Two-column table of intra-protein aligning peptide pairs.
#' @param verbose Report progress every 2,000 proteins.
#' @return A `data.table` with `pro_id` followed by one column per sample.
#' @export
polyclonal_scores <- function(hits, pro_id, pairs, verbose = FALSE) {
  hits <- data.table::as.data.table(hits)
  if (length(pro_id) != nrow(hits))
    stop(sprintf("pro_id has length %d but hits has %d rows.",
                 length(pro_id), nrow(hits)), call. = FALSE)

  sample_cols <- names(hits)[-1L]
  peptides <- as.character(hits[[1L]])
  mapped <- !is.na(pro_id)
  proteins <- unique(pro_id[mapped])

  if (length(proteins) == 0L) {
    out <- data.table::data.table(pro_id = character())
    for (j in sample_cols) data.table::set(out, j = j, value = numeric())
    return(out[])
  }

  # Alignment pairs restricted to peptides actually present, as an index pair
  # list so no string matching happens inside the per-protein loop.
  pairs <- data.table::as.data.table(pairs)
  if (ncol(pairs) >= 2L && nrow(pairs) > 0L) {
    a <- match(as.character(pairs[[1L]]), peptides)
    b <- match(as.character(pairs[[2L]]), peptides)
    ok <- !is.na(a) & !is.na(b) & a != b
    a <- a[ok]; b <- b[ok]
  } else {
    a <- integer(0); b <- integer(0)
  }

  rows_by_protein <- split(which(mapped), factor(pro_id[mapped], levels = proteins))
  # Which alignment edges belong to which protein (both endpoints share one).
  edge_protein <- if (length(a)) pro_id[a] else character(0)
  edges_by_protein <- if (length(a))
    split(seq_along(a), factor(edge_protein, levels = proteins)) else
      stats::setNames(vector("list", length(proteins)), proteins)

  hit_mat <- as.matrix(hits[, sample_cols, with = FALSE])
  scores <- matrix(0, nrow = length(proteins), ncol = length(sample_cols),
                   dimnames = list(NULL, sample_cols))

  for (p in seq_along(proteins)) {
    if (verbose && p %% 2000L == 0L)
      message(sprintf("  polyclonal: %d / %d proteins", p, length(proteins)))

    rows <- rows_by_protein[[p]]
    if (length(rows) == 0L) next
    block <- hit_mat[rows, , drop = FALSE]

    e <- edges_by_protein[[p]]
    # Fast path: no internal alignments, or a single peptide.
    if (is.null(e) || length(e) == 0L || length(rows) == 1L) {
      scores[p, ] <- colSums(block, na.rm = TRUE)
      next
    }

    # Re-express this protein's edges as local row indices 1..length(rows).
    ea <- match(a[e], rows); eb <- match(b[e], rows)
    keep <- !is.na(ea) & !is.na(eb)
    if (!any(keep)) {
      scores[p, ] <- colSums(block, na.rm = TRUE)
      next
    }
    g <- igraph::graph_from_edgelist(cbind(ea[keep], eb[keep]), directed = FALSE)
    if (igraph::vcount(g) < length(rows))
      g <- igraph::add_vertices(g, length(rows) - igraph::vcount(g))

    memo <- new.env(parent = emptyenv())
    for (s in seq_along(sample_cols)) {
      hv <- which(block[, s] == 1)
      if (length(hv) == 0L) next
      if (length(hv) == 1L) { scores[p, s] <- 1; next }

      key <- paste(hv, collapse = ",")
      cached <- memo[[key]]
      if (!is.null(cached)) { scores[p, s] <- cached; next }

      sub <- igraph::induced_subgraph(g, hv)
      val <- if (igraph::ecount(sub) == 0L) length(hv) else .independent_hits(sub)
      memo[[key]] <- val
      scores[p, s] <- val
    }
  }

  out <- data.table::data.table(pro_id = proteins)
  for (s in seq_along(sample_cols))
    data.table::set(out, j = sample_cols[s], value = scores[, s])
  out[]
}
