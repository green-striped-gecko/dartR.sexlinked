# Slow file: skipped on CRAN to keep the check time under 10 minutes; runs
# with devtools::test() and on GitHub Actions.
skip_on_cran()

# File-backed objects: the per-locus counts by sex (gl.report.sexlinked,
# gl.filter.sexlinked, gl.keep.sexlinked) and the per-individual counts of
# gl.infer.sex are read in blocks of loci (dartR.base::utils.fbm.colwise,
# utils.fbm.rowwise, utils.fbm.reader). Results and messages must be
# identical to those of the gen-backed object for any block size, and no
# call may decode the whole FBM.

lbp <- LBP[, 1:600]
lbp.f <- suppressWarnings(gl.gen2fbm(lbp, verbose = 0))

# genlights (also inside lists) are compared by their decoded genotypes and
# metadata, as a file-backed input gives file-backed outputs; plots dropped
strip <- function(v) {
  if (methods::is(v, "genlight")) {
    return(list(m = as.matrix(v), ind = indNames(v), loc = locNames(v),
                pop = pop(v), lm = v@other$loc.metrics,
                im = v@other$ind.metrics))
  }
  if (is.list(v) && !is.data.frame(v)) {
    v[vapply(v, inherits, logical(1), "ggplot")] <- NULL
    v <- lapply(v, strip)
  }
  v
}

# value, console output and messages of a call
run <- function(FUN, x) {
  msg <- character(0)
  out <- utils::capture.output(val <- withCallingHandlers(
    suppressWarnings(FUN(x)),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  list(val = strip(val), out = c(out, msg))
}

calls <- list(
  report = function(x) gl.report.sexlinked(x, system = "xy",
                                           plot.display = FALSE, verbose = 3),
  report.zw = function(x) gl.report.sexlinked(x, system = "zw",
                                              plot.display = FALSE,
                                              verbose = 3),
  filter = function(x) gl.filter.sexlinked(x, system = "xy",
                                           plot.display = FALSE, verbose = 3),
  keep = function(x) gl.keep.sexlinked(x, system = "xy",
                                       plot.display = FALSE, verbose = 3)
)

# all three sex assignments (W/Y, Z/X and gametologs) need the full data
infer <- function(x, system) {
  k <- suppressWarnings(gl.keep.sexlinked(x, system = system,
                                          plot.display = FALSE, verbose = 0))
  gl.infer.sex(k, system = system, seed = 100, verbose = 3)
}

test_that("gen- and file-backed objects give identical results", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  for (bs in c(1, 1e-6)) {                   # 1e-6: one locus per block
    options(bigstatsr.block.sizeGB = bs)
    for (nm in names(calls)) {
      g <- run(calls[[nm]], lbp)
      f <- run(calls[[nm]], lbp.f)
      expect_identical(f$val, g$val, label = paste(nm, bs))
      expect_identical(f$out, g$out, label = paste(nm, bs, "output"))
    }
  }
})

test_that("gl.infer.sex gives identical results for a file-backed object", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  f.lbp <- suppressWarnings(gl.gen2fbm(LBP, verbose = 0))
  for (system in c("xy", "zw")) {
    g <- run(function(x) infer(x, system), LBP)
    f <- run(function(x) infer(x, system), f.lbp)
    expect_identical(f$val, g$val, label = system)
    expect_identical(f$out, g$out, label = paste(system, "output"))
  }
  # with zw the gametolog assignment ran (it needs at least 5 gametologs)
  expect_false(all(is.na(g$val$gametolog.sex)))
})

test_that("file-backed objects are never decoded as a whole", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  }, .package = "dartR.base")
  for (nm in c("report", "report.zw")) {
    expect_no_error(suppressMessages(suppressWarnings(
      utils::capture.output(calls[[nm]](lbp.f)))), message = nm)
  }
})
