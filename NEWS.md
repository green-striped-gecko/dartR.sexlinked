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
