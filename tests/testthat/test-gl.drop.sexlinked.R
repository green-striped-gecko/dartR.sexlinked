# Characterization test for gl.drop.sexlinked() (deprecated in favour of
# gl.filter.sexlinked()). Captured on commit dc41e43 (version 1.2.6), then
# updated only for the changes approved in function-review/reports/
# dartR.sexlinked/gl.drop.sexlinked.md. A failing expectation here means
# behaviour changed; it does not mean the old behaviour was correct.

test_that("deprecation warning names the replacement", {
  expect_warning(gl.drop.sexlinked(LBP, system = "xy", plot.display = FALSE,
                                   verbose = 0),
                 "gl.filter.sexlinked")
})

test_that("returns what gl.filter.sexlinked returns", {
  # change 1: previously returned the warning text (a character string)
  r <- suppressWarnings(gl.drop.sexlinked(LBP, system = "xy",
                                          plot.display = FALSE, verbose = 0))
  f <- gl.filter.sexlinked(LBP, system = "xy", plot.display = FALSE,
                           verbose = 0)
  expect_s4_class(r, "genlight")
  expect_equal(nLoc(r), nLoc(f))
  expect_identical(locNames(r), locNames(f))
  expect_identical(as.matrix(r), as.matrix(f))
})

test_that("pre-deprecation arguments are accepted", {
  # change 1: previously "unused arguments"
  r <- suppressWarnings(gl.drop.sexlinked(LBP, system = "xy", ncores = 1,
                                          plot.display = FALSE, verbose = 0))
  expect_s4_class(r, "genlight")
})

test_that("history records the caller's arguments", {
  # addendum A1: previously recorded x = x, system = system
  r <- suppressWarnings(gl.drop.sexlinked(LBP, system = "xy",
                                          plot.display = FALSE, verbose = 0))
  f <- gl.filter.sexlinked(LBP, system = "xy", plot.display = FALSE,
                           verbose = 0)
  expect_identical(tail(r@other$history, 1), tail(f@other$history, 1))
})
