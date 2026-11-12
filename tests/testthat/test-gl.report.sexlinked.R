# Characterization test for gl.report.sexlinked().
# Captured on commit f293fba (version 1.2.2), then updated only for the
# changes approved in function-review/reports/dartR.sexlinked/
# gl.report.sexlinked.md (changes 1-6). A failing expectation here means
# behaviour changed; it does not mean the old behaviour was correct.

run_report <- function(x, system = "xy", ...) {
  gl.report.sexlinked(x, system = system, plot.display = FALSE,
                      verbose = 0, ...)
}

count_categories <- function(r, system) {
  first <- if (system == "xy") "y.linked" else "w.linked"
  third <- if (system == "xy") "x.linked" else "z.linked"
  c(first = sum(r[[first]]), sex.biased = sum(r$sex.biased),
    third = sum(r[[third]]), gametolog = sum(r$gametolog))
}

test_that("LBP, xy: table shape, categories and key values", {
  x0 <- LBP
  r <- run_report(LBP, "xy")
  expect_identical(LBP, x0) # input object untouched
  expect_s3_class(r, "data.frame")
  expect_equal(dim(r), c(1000L, 23L))
  expect_equal(names(r)[11], "y.linked")
  expect_equal(unname(count_categories(r, "xy")), c(1, 9, 66, 1))
  expect_equal(sum(r$p.adjusted <= 0.01), 10L)
  expect_equal(sum(r$count.F.scored + r$count.F.miss == 162), 1000L)
  expect_equal(sum(r$count.M.scored + r$count.M.miss == 211), 1000L)
  expect_equal(sum(is.na(r$stat.p.value)), 10L)
})

test_that("LBP, zw: categories", {
  r <- run_report(LBP, "zw")
  expect_equal(names(r)[11], "w.linked")
  expect_equal(unname(count_categories(r, "zw"))[c(1, 3)], c(0, 1))
})

test_that("testset.gl and platypus.gl: categories", {
  r <- run_report(testset.gl, "xy")
  expect_equal(dim(r), c(255L, 23L))
  expect_equal(unname(count_categories(r, "xy")), c(0, 0, 0, 0))
  # platypus.gl has both 'Sex' and 'sex' columns; the first match is used
  r <- suppressWarnings(run_report(platypus.gl, "xy"))
  expect_equal(unname(count_categories(r, "xy")), c(1, 4, 1, 2))
})

test_that("serial and parallel results agree", {
  skip_on_cran()
  r1 <- run_report(LBP, "xy")
  r2 <- run_report(LBP, "xy", ncores = 2)
  expect_equal(r1, r2, ignore_attr = TRUE)
})

test_that("Fisher's test uses the observed counts (change 6)", {
  r <- run_report(LBP, "xy")
  i <- which(r$count.F.miss == 0 & r$count.M.miss > 0)[1]
  obs <- matrix(unlist(r[i, c("count.F.miss", "count.M.miss",
                              "count.F.scored", "count.M.scored")]), 2)
  expect_equal(r$p.value[i], fisher.test(obs)$p.value)
  expect_equal(r$ratio[i], unname(fisher.test(obs)$estimate))
})

test_that("edge cases", {
  x <- testset.gl
  full <- run_report(x, "xy")
  # Exactly one female runs (change 3)
  x1 <- x
  f <- which(x1@other$ind.metrics$sex == "Female")
  x1@other$ind.metrics$sex[f[-1]] <- "Unknown"
  r <- run_report(x1, "xy")
  expect_true(all(r$count.F.scored + r$count.F.miss == 1))
  # No males errors (change 2)
  x2 <- x
  x2@other$ind.metrics$sex[x2@other$ind.metrics$sex == "Male"] <- "Unknown"
  expect_error(run_report(x2, "xy"), "114 females and 0 males")
  # Sex is read by row, so a missing or renamed 'id' column changes
  # nothing (change 1)
  x3 <- x
  x3@other$ind.metrics$id <- NULL
  expect_equal(run_report(x3, "xy"), full)
  x4 <- x
  x4@other$ind.metrics$id <- paste0("s", seq_len(nInd(x4)))
  expect_equal(run_report(x4, "xy"), full)
  # SilicoDArT is rejected (change 5)
  expect_error(run_report(testset.gs, "xy"), "found SilicoDArT")
  # Missing system
  expect_error(run_report(x, NULL), "sex-determination system")
})

test_that("parallel workers are stopped after an error (change 4)", {
  skip_on_cran()
  n <- nrow(showConnections())
  x2 <- testset.gl
  x2@other$ind.metrics$sex[x2@other$ind.metrics$sex == "Male"] <- "Unknown"
  expect_error(run_report(x2, "xy", ncores = 2))
  expect_equal(nrow(showConnections()), n)
})
