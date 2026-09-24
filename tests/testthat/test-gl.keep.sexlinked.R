# Characterization test for gl.keep.sexlinked().
# Captured on commit f293fba (version 1.2.2), then updated only for the
# changes approved in function-review/reports/dartR.sexlinked/
# gl.keep.sexlinked.md. A failing expectation here means
# behaviour changed; it does not mean the old behaviour was correct.

run_keep <- function(x, system = "xy", ...) {
  gl.keep.sexlinked(x, system = system, plot.display = FALSE,
                    verbose = 0, ...)
}

n_loc <- function(o) if (is.null(o)) 0L else nLoc(o)

test_that("LBP, xy: list structure, subsets and metadata", {
  x0 <- LBP
  r <- run_keep(LBP, "xy")
  expect_identical(LBP, x0) # input object untouched
  expect_named(r, c("results.table", "y.linked", "sex.biased", "x.linked",
                    "gametolog"))
  expect_equal(dim(r$results.table), c(1000L, 23L))
  expect_equal(unname(sapply(r[-1], n_loc)), c(1L, 9L, 66L, 1L))
  for (g in r[-1]) {
    expect_equal(nInd(g), 376L)
    expect_equal(unique(ploidy(g)), 2L)
    expect_equal(nrow(g@other$loc.metrics), nLoc(g))
  }
  idx <- r$results.table$index[r$results.table$x.linked]
  expect_identical(locNames(r$x.linked), locNames(LBP)[idx])
  # One history entry is added to each returned object (change 13)
  expect_equal(length(r$x.linked@other$history),
               length(LBP@other$history) + 1L)
})

test_that("empty categories return NULL", {
  r <- run_keep(LBP, "zw")
  expect_named(r, c("results.table", "w.linked", "sex.biased", "z.linked",
                    "gametolog"))
  expect_null(r$w.linked)
  expect_equal(unname(sapply(r[-1], n_loc)), c(0L, 10L, 1L, 66L))
  r <- suppressWarnings(run_keep(platypus.gl, "xy"))
  expect_equal(unname(sapply(r[-1], n_loc)), c(1L, 4L, 1L, 2L))
})

test_that("results table equals gl.report.sexlinked()", {
  expect_equal(run_keep(LBP, "xy")$results.table,
               gl.report.sexlinked(LBP, "xy", plot.display = FALSE,
                                   verbose = 0))
})

test_that("serial and parallel results agree", {
  skip_on_cran()
  # The history entries differ (they record ncores), so they are removed
  no_history <- function(r) {
    lapply(r, function(o) {
      if (inherits(o, "genlight")) o@other$history <- NULL
      o
    })
  }
  expect_equal(no_history(run_keep(LBP, "xy")),
               no_history(run_keep(LBP, "xy", ncores = 2)))
})

test_that("edge cases", {
  x <- testset.gl
  # Plain genlight: loc.metrics follow the kept loci (change 12)
  r <- run_keep(as(LBP, "genlight"), "xy")
  expect_equal(nLoc(r$x.linked), 66L)
  expect_equal(nrow(r$x.linked@other$loc.metrics), 66L)
  idx <- r$results.table$index[r$results.table$x.linked]
  expect_identical(r$x.linked@other$loc.metrics$AlleleID,
                   LBP@other$loc.metrics$AlleleID[idx])
  # Exactly one female runs (change 3)
  x1 <- x
  f <- which(x1@other$ind.metrics$sex == "Female")
  x1@other$ind.metrics$sex[f[-1]] <- "Unknown"
  expect_type(run_keep(x1, "xy"), "list")
  # No males errors (change 2)
  x2 <- x
  x2@other$ind.metrics$sex[x2@other$ind.metrics$sex == "Male"] <- "Unknown"
  expect_error(run_keep(x2, "xy"), "114 females and 0 males")
  # Sex is read by row, so a missing 'id' column changes nothing (change 1)
  x3 <- x
  x3@other$ind.metrics$id <- NULL
  expect_equal(run_keep(x3, "xy")$results.table,
               run_keep(x, "xy")$results.table)
  # SilicoDArT is rejected (change 5)
  expect_error(run_keep(testset.gs, "xy"), "found SilicoDArT")
  expect_error(run_keep(x, NULL), "sex-determination system")
})
