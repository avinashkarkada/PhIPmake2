# phipmake2

[![R-CMD-check](https://github.com/avinashkarkada/PhIPmake2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/avinashkarkada/PhIPmake2/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](LICENSE)

Post-alignment summarisation for [PhIP-Seq](https://www.nature.com/articles/s41596-018-0025-6)
data. Takes peptide-by-sample matrices and produces per-library splits,
annotated tables, and protein-level collapses.

A drop-in replacement for the `phipmake` step of the Larman Lab pipeline: same
inputs, same output names, same format.

## What it does

Given `counts.csv`, `fold_change.csv`, `enrichment.csv` and `Hits.csv` for a
screen:

- splits each matrix into its constituent peptide sub-libraries
- joins peptide and protein annotations
- collapses peptides to proteins by sum (`prosum`), maximum (`promax`), and
  independent-epitope count (`polyclonal`)
- writes hit-masked variants of counts, fold change and enrichment

Outputs are written per sub-library and as whole-screen files, annotated and
unannotated.

## Installation

```r
# install.packages("remotes")
remotes::install_github("avinashkarkada/PhIPmake2")
```

Needs R 4.0 or newer, `data.table` 1.14 or newer and `igraph` 1.2 or newer.

On a cluster, point `R_LIBS_USER` at a writable directory of your own followed
by whichever shared library provides `data.table` and `igraph`:

```bash
module load r/4.3.0
export R_LIBS_USER=/your/lib:/path/to/shared/R/library
Rscript -e 'remotes::install_github("avinashkarkada/PhIPmake2")'
```

## Usage

```r
library(phipmake2)

run_phipmake2(
  wd            = "/path/to/screen",
  stages        = c("counts", "foldchange", "enrichment", "hits"),
  metadata_path = "/path/to/PeptideLibraries"
)
```

A `drake_params.tsv` already in the screen directory is picked up
automatically, and dash-separated plan strings work:

```r
run_phipmake2(wd, stages = "Counts-FoldChange-Enrichment-Hits-Polyclonal")
```

The metadata directory can also come from the `PHIPMAKE2_METADATA` environment
variable rather than being passed each time.

### Command line

```bash
Rscript "$(Rscript -e 'cat(system.file("cli/phipmake2.R", package="phipmake2"))')" \
  --wd /path/to/screen \
  --stages Counts-FoldChange-Enrichment-Hits \
  --metadata /path/to/PeptideLibraries
```

### SLURM

```bash
sbatch --export=wd="/path/to/screen/",plan="Counts-FoldChange-Enrichment-Hits" \
  inst/slurm/runphipmake2.sh
```

## Performance

On a full screen (Counts-FoldChange-Enrichment-Hits, 2 CPUs each), PhIPmake2 ran
in 2m 39s against phipmake's 4h 12m, with peak resident memory of 3.9 GB against
38.9 GB. All 112 peptide-level outputs were byte-identical.

Collapsing peptides to proteins is the expensive step. Done as a loop over
proteins it rescans every peptide once per protein and copies the output frame
on each iteration; phipmake2 does it as a single grouped `data.table` pass.

On a synthetic matrix at the scale of a large sub-library (106,679 peptides,
145 samples), with identical output either way:

| | time |
|---|---|
| loop over proteins | 126.58 s |
| grouped pass | 1.04 s |

`bench/bench_collapse.R` runs the comparison; pass `--full` to include the
largest case.

Peak memory is bounded by processing one input matrix and one sub-library at a
time, and by appending whole-screen files to disk as each sub-library finishes
rather than assembling them in memory.

## Peptides with no annotation

A peptide absent from its library's annotation file has no protein, so it
cannot contribute to a protein-level total. phipmake2 excludes those peptides
and reports how many there were. Their peptide-level data is untouched.

This matters because the peptide-to-protein map has to stay the same length as
the data it describes. Dropping the unannotated entries first produces a shorter
vector, and R recycles it when it is used for grouping, quietly shifting
peptides into neighbouring proteins. One missing peptide shifts everything after
it; a library missing much of its annotation ends up with totals bearing no
relationship to the protein they are filed under. `collapse_peptides()` rejects
a mis-sized map rather than accepting the recycle.

Excluding those peptides makes a protein-level result *incomplete*, not wrong.
Completeness needs the annotation file regenerated, which is a metadata job.

## Comparing against a previous run

```r
compare_outputs(old = "previous_run", new = "new_run")
```

Values and bytes are reported separately, so "same answer, written differently"
is distinguishable from "wrong answer".

`inst/validate/validate_screen.R` wraps this into one job: it checks the
metadata, runs the pipeline into a fresh directory with timing and peak memory,
then compares against the existing output.

See [MIGRATION.md](MIGRATION.md) for what to expect when switching a screen
over.

## Scope

This covers the summarisation step only. Alignment, library splitting, EdgeR
enrichment, hit calling, ARscore and ARscape are separate and unchanged.

AVARDA and epitopefindr are accepted in a plan string but not implemented; they
are skipped with a warning. AVARDA is planned.

## Development

```r
devtools::test()   # 143 tests
devtools::check()
```

## Acknowledgements

Supersedes `phipmake` by Brandon Sie, whose output format and analysis
definitions this follows. The polyclonal independence filter is Daniel Monaco's.
