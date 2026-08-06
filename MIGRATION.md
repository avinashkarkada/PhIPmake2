# Switching a screen over from phipmake

## What stays the same

**Inputs.** Same four matrices in the same working directory (`counts.csv`,
`fold_change.csv`, `enrichment.csv`, `Hits.csv`) and the same
`drake_params.tsv`. Nothing upstream changes.

**Outputs.** Same file names, same locations, same columns:

| | pattern |
|---|---|
| whole screen | `{screen}_{stem}.tsv`, `{screen}_{stem}_annotated.tsv` |
| per library | `{lib}/{screen}_{lib}_{stem}.tsv`, `{lib}/{screen}_{lib}_{stem}_annotated.tsv` |

`stem` is one of `Counts`, `Counts_prosum`, `FoldChange`, `EdgeR`,
`EdgeR_promax`, `Hits`, `Hits_prosum`, `Hits_counts`, `Hits_foldchange`,
`Hits_foldchange_promax`, `Hits_enrichment`, `Hits_enrichment_promax`,
`Polyclonal`. The `Counts`, `FoldChange` and `EdgeR` parts follow
`counts_type`, `foldchange_type` and `enrichment_type` from `drake_params.tsv`.

Tab separated, `NA` for missing, character columns quoted. Column order is
identifier, annotation fields, samples. Peptides keep input order within a
sub-library; proteins appear in order of first appearance.

**Downstream.** ARscore matches `^(.*)_FoldChange_annotated(.*)$` and takes
`_Hits_foldchange_annotated.tsv`. ARscape takes the same two. Both are
peptide-level and unchanged.

## What you need to change

Only the job submission.

1. Point the dispatcher at `runphipmake2.sh` instead of `runphipmake.sh`
2. Drop the `clean=` and `njobs=` exports; add `metadata=` and `threads=`
3. Load an R 4.x module instead of the R 3.6 one, with `R_LIBS_USER` pointing
   at a library built for it
4. Install phipmake2 into that library

phipmake2 requires R 4.0 or newer, so it cannot live in an R 3.6 library
alongside phipmake. The two packages have different names and no overlapping
function names, so keeping both around causes no conflict.

## What will differ in the output

**Protein-level numbers can change.** `prosum`, `promax` and `Polyclonal` are
computed from the peptide-to-protein map. Previously that map had its `NA`s
dropped before use, which shortened it relative to the data; R then recycled it
during grouping, shifting peptides into neighbouring proteins. If a library's
annotation covered every peptide, nothing changes. Otherwise the new numbers
are correct and the old ones were not.

Peptide-level files are unaffected.

**Some conditions now stop the run instead of continuing.** A library with no
annotation directory, an ambiguous annotation file, or a peptide identifier that
does not parse are all errors now. Previously these resolved silently to
whichever file matched first alphabetically, or produced a malformed identifier
table.

**Empty protein groups produce `NA`** rather than `-Inf`.

**Unannotated peptides are counted** in a warning at the end of each run.

**AVARDA, epitopefindr and pairwise are skipped** with a warning if they appear
in a plan string, rather than aborting the run.

**There is no `.drake` cache.** `drake::clean()` and the `clean=` variable have
no equivalent and are not needed, since rerunning is cheap.

**`njobs` is gone.** Parallelism is `data.table`'s internal threading, set with
`threads=` or `--threads`.

## Suggested procedure

1. Run phipmake2 into a fresh directory, not over existing output
2. `compare_outputs(old, new)`
3. Confirm peptide-level files come back `identical`
4. For any `values_differ` on a protein-level file, check the run log for the
   unannotated-peptide warning, which is the expected cause
5. Confirm ARscore and ARscape still run
