# dartR.sexlinked 1.2.6

## gl.drop.sexlinked()

* The deprecated `gl.drop.sexlinked()` now warns and then runs
  `gl.filter.sexlinked()` with the same arguments, returning the autosomal
  genlight. Since March 2026 it returned its warning text, so
  `x <- gl.drop.sexlinked(x, "xy")` replaced the data with a character
  string, and calls with `ncores`, `plot.display` or `verbose` failed with
  "unused arguments".

## gl.report.sexlinked()

* Sex is now read from the rows of `ind.metrics` in individual order instead
  of matching an `id` column against the individual names. Objects whose `id`
  column was missing or differed from `indNames()` previously got all-zero
  counts and every locus reported as autosomal; they now get real results.
* The function now stops with an error when either sex has no individuals,
  instead of reporting every locus as autosomal.
* A sex with a single individual no longer causes an error.
* SilicoDArT (presence/absence) data are no longer accepted, because
  heterozygosity is undefined for them.
* Fisher's exact tests now use the observed counts. Previously every zero
  count was replaced by 1, which inflated p-values when few individuals were
  sexed. P-values, `ratio` and `stat` change, and some loci change category
  (on a 200-individual, 10,000-SNP platypus dataset: 85 to 95 x-linked loci).
  The chi-square branch (1,000 or more individuals) is unchanged.
* `plot.theme` is now applied to both plots.
* Parallel workers are stopped when the function exits with an error.
* Classification loops are vectorised: about 2 times faster in serial runs.

## gl.keep.sexlinked()

* Receives the classification fixes made to `gl.report.sexlinked()` above:
  sex read by row, an error when a sex has no individuals, no crash with a
  single individual of a sex, SNP data only, Fisher's test on the observed
  counts, `plot.theme` applied, workers stopped on error, vectorised loops.
  Its `results.table` is identical to the output of `gl.report.sexlinked()`.
* Returned genlight objects now carry only their own loci's `loc.metrics`
  when the input is a plain genlight (previously all loci's metrics).
* Each returned genlight object records the call in `@other$history`.
* Empty categories are documented as `NULL`.

## gl.filter.sexlinked()

* Receives the same classification fixes as `gl.report.sexlinked()` and
  `gl.keep.sexlinked()`. Datasets with no males or an unmatched `id` column
  previously came back unfiltered; they now error or are filtered correctly.
  The loci returned are exactly those not returned by `gl.keep.sexlinked()`.
* The returned object carries only its own loci's `loc.metrics` when the
  input is a plain genlight, and records the call in `@other$history`.
* Returning NULL when every locus is sex-linked is now documented.

## gl.infer.sex()

* `agreed.sex` changes in two cases. An individual with no usable genotypes
  is now `NA` (was `*F`). A tie between one `F` and one `M` assignment is now
  `*?` (was always `*F`, whichever way round the tie was; 28 of 376
  individuals in LBP).
* The function now stops with an error when `gl_sexlinked` is not the output
  of `gl.keep.sexlinked()`, was built with the other `system`, or contains no
  sex-linked loci. Previously a mismatched `system` returned `*F` for every
  individual with ids replaced by row numbers.
* New `verbose` argument; messages follow the standard dartR verbosity levels.
* The caller's random number stream is restored on exit. Previously the
  `seed` used for k-means replaced it.
* Individuals are matched by position, so duplicated individual names no
  longer set every X-/Z-linked and gametolog assignment to `NA`.
* Documentation now states that only the top 5 gametologs are used, the
  minimum numbers of loci, and that k-means assumes both sexes are present.
