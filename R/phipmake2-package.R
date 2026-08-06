#' phipmake2: fast, memory-bounded PhIP-Seq summarisation
#'
#' A rewrite of the Larman Lab phipmake post-alignment stage. See
#' [run_phipmake2()] for the entry point, and the package README for the
#' rationale behind the rewrite.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats setNames
#' @importFrom utils head packageVersion
"_PACKAGE"

# Columns and symbols referenced by data.table non-standard evaluation.
utils::globalVariables(c("n", "base", "version", "pro_id", "verdict", "bytes", ".SD"))
