# dartR.sexlinked 1.2.2.9000

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
