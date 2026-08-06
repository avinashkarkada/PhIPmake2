# phipmake2 2.0.0

First release. Supersedes `phipmake` 0.2.8, keeping its output format and
analysis definitions. See [MIGRATION.md](MIGRATION.md) for what changes when
switching an existing screen over.

## Analysis

* The peptide-to-protein map is kept at full length, with `NA` for peptides that
  have no annotation, and those rows are excluded explicitly. Dropping the `NA`s
  first shortens the map relative to the data, which R recycles during grouping
  and which silently assigns peptides to the wrong proteins. This affects
  `prosum`, `promax` and `Polyclonal`.
* Unannotated peptides are counted and reported in a warning.
* Empty protein groups yield `NA` rather than `-Inf`.

## Robustness

* Library names are matched exactly rather than by substring, so one library
  name is no longer also matched by a longer name sharing its prefix.
* Peptide identifiers are parsed with an anchored regex that tolerates
  underscores in both the library base name and the sequence.
* Annotation files are resolved by exact library name; a missing or ambiguous
  match is an error.
* Hit masking checks that the data and hit matrices share a row order before
  masking by position.

## Performance

* Protein collapse is a single grouped `data.table` pass instead of a loop over
  proteins, with the same output. See `bench/bench_collapse.R`.
* Polyclonal scoring builds each protein's alignment graph once and induces it
  per sample, memoises repeated hit sets, and skips the graph entirely for
  proteins with no internal alignments.
* Peptide identifiers are parsed once per run rather than once per data file.

## Memory

* One input matrix and one sub-library are held at a time.
* Whole-screen files are appended to disk per sub-library rather than assembled
  in memory.
* Hit masking mutates column by column rather than copying the matrix.
* No cache on disk.

## Interface

* `run_phipmake2()` is the entry point; there is no plan to build first.
* `compare_outputs()` diffs two output directories.
* CLI at `inst/cli/phipmake2.R`, SLURM script at `inst/slurm/runphipmake2.sh`.
* AVARDA, epitopefindr and pairwise are accepted in a plan string and skipped
  with a warning. AVARDA is planned.
* Depends only on `data.table` and `igraph`. Requires R 4.0 or newer.
* 126 tests.
