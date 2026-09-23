# Characterization test for gl.filter.sexlinked().
# Captured on commit f293fba (version 1.2.2), then updated only for the
# changes approved in function-review/reports/dartR.sexlinked/
# gl.filter.sexlinked.md (dartR.base). A failing expectation here means
# behaviour changed; it does not mean the old behaviour was correct.

run_filter <- function(x, system = "xy", ...) {
  gl.filter.sexlinked(x, system = system, plot.display = FALSE,
                      verbose = 0, ...)
}

autosomal_index <- function(x, system) {
  t <- gl.report.sexlinked(x, system = system, plot.display = FALSE,
                           verbose = 0)
  t$index[rowSums(t[, sapply(t, is.logical)]) == 0]
}

test_that("LBP: returns the autosomal loci with metadata", {
  x0 <- LBP
  r <- run_filter(LBP, "xy")
  expect_identical(LBP, x0) # input object untouched
  expect_s4_class(r, "genlight")
  expect_equal(nInd(r), 376L)
  expect_equal(nLoc(r), 923L)
  expect_equal(unique(ploidy(r)), 2L)
  expect_equal(nrow(r@other$loc.metrics), 923L)
  expect_identical(locNames(r), locNames(LBP)[autosomal_index(LBP, "xy")])
  # One history entry is added (change 13)
  expect_equal(length(r@other$history), length(LBP@other$history) + 1L)
  expect_equal(nLoc(run_filter(LBP, "zw")), 923L)
})

test_that("platypus.gl and testset.gl", {
  expect_equal(nLoc(suppressWarnings(run_filter(platypus.gl, "xy"))), 992L)
  expect_equal(nLoc(run_filter(testset.gl, "xy")), 255L)
})

test_that("serial and parallel results agree", {
  skip_on_cran()
  # The history entries differ (they record ncores), so they are removed
  no_history <- function(o) {
    o@other$history <- NULL
    o
  }
  expect_equal(no_history(run_filter(LBP, "xy")),
               no_history(run_filter(LBP, "xy", ncores = 2)))
})

test_that("all loci sex-linked returns NULL", {
  aut <- autosomal_index(LBP, "xy")
  sl <- gl.keep.loc(LBP, locNames(LBP)[-aut], verbose = 0)
  expect_null(run_filter(sl, "xy"))
})

test_that("filter and keep partition the loci", {
  k <- gl.keep.sexlinked(LBP, "xy", plot.display = FALSE, verbose = 0)
  sl <- unlist(lapply(k[-1], function(o) if (is.null(o)) NULL else locNames(o)))
  expect_identical(locNames(run_filter(LBP, "xy")),
                   setdiff(locNames(LBP), sl))
})

test_that("edge cases", {
  x <- testset.gl
  # Plain genlight: loc.metrics follow the kept loci (change 12)
  r <- run_filter(as(LBP, "genlight"), "xy")
  expect_equal(nLoc(r), 923L)
  expect_equal(nrow(r@other$loc.metrics), 923L)
  # Exactly one female runs (change 3)
  x1 <- x
  f <- which(x1@other$ind.metrics$sex == "Female")
  x1@other$ind.metrics$sex[f[-1]] <- "Unknown"
  expect_s4_class(run_filter(x1, "xy"), "genlight")
  # No males errors instead of returning every locus (change 2)
  x2 <- x
  x2@other$ind.metrics$sex[x2@other$ind.metrics$sex == "Male"] <- "Unknown"
  expect_error(run_filter(x2, "xy"), "114 females and 0 males")
  # Sex is read by row, so a missing 'id' column changes nothing (change 1)
  x3 <- x
  x3@other$ind.metrics$id <- NULL
  expect_identical(locNames(run_filter(x3, "xy")),
                   locNames(run_filter(x, "xy")))
  # SilicoDArT is rejected (change 5)
  expect_error(run_filter(testset.gs, "xy"), "found SilicoDArT")
  expect_error(run_filter(x, NULL), "sex-determination system")
})
