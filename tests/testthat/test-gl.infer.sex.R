# Characterization test for gl.infer.sex().
# Captured on commit 179080b (version 1.2.2) before review, then updated only
# for the changes approved in function-review/reports/dartR.sexlinked/
# gl.infer.sex.md (dartR.base). A failing expectation here means behaviour
# changed; it does not mean the old behaviour was correct.

keep_quiet <- function(x, system) {
  suppressWarnings(gl.keep.sexlinked(x, system = system,
                                     plot.display = FALSE, verbose = 0))
}
infer_quiet <- function(...) gl.infer.sex(..., verbose = 0)

k_xy <- keep_quiet(LBP, "xy")
k_zw <- keep_quiet(LBP, "zw")

test_that("LBP, xy: output shape and assignments", {
  r <- infer_quiet(k_xy, system = "xy", seed = 100)
  expect_s3_class(r, "data.frame")
  expect_equal(dim(r), c(376L, 11L))
  expect_named(r, c("id", "y.linked.sex", "#called", "#missing",
                    "x.linked.sex", "#Het.x", "#Hom.x", "gametolog.sex",
                    "#Het.g", "#Hom.g", "agreed.sex"))
  expect_identical(r$id, indNames(LBP))
  # The 28 one-against-one ties were '*F' before change 3
  expect_equal(as.vector(table(r$agreed.sex)), c(28L, 136L, 212L))
  expect_equal(names(table(r$agreed.sex)), c("*?", "F", "M"))
  expect_identical(row.names(r), indNames(LBP))
  expect_equal(unname(unlist(r[2, 2:7])), c("M", "1", "0", "M", "2", "56"))
  # One Y-linked locus, one gametolog (< 5): gametolog columns are NA
  expect_true(all(is.na(r$gametolog.sex)))
})

test_that("LBP, zw: gametologs only", {
  r <- infer_quiet(k_zw, system = "zw", seed = 100)
  expect_true(all(is.na(r$w.linked.sex)))
  expect_true(all(is.na(r$z.linked.sex)))
  expect_equal(as.vector(table(r$agreed.sex)), c(101L, 275L))
  # Only the top 5 gametologs are counted
  expect_true(all(r$`#Het.g` + r$`#Hom.g` <= 5))
})

test_that("input is checked (change 1)", {
  # Before: all individuals '*F' and ids replaced by row numbers
  expect_error(infer_quiet(k_xy, system = "zw", seed = 100),
               "built with system = 'xy'")
  expect_error(infer_quiet(k_zw, system = "xy", seed = 100),
               "built with system = 'zw'")
  # Before: "invalid 'times' argument"
  k0 <- k_xy
  for (s in c("y.linked", "sex.biased", "x.linked", "gametolog")) {
    k0[s] <- list(NULL)
  }
  expect_error(infer_quiet(k0, system = "xy", seed = 1),
               "No sex-linked loci")
  expect_error(infer_quiet(LBP, system = "xy", seed = 1),
               "must be the list returned by gl.keep.sexlinked")
})

test_that("no information gives NA (change 2)", {
  gx <- k_xy$x.linked
  m <- as.matrix(gx)
  m[1, ] <- NA
  kna <- k_xy
  kna$x.linked <- new("genlight", m, ploidy = 2L, ind.names = indNames(gx),
                      loc.names = locNames(gx))
  kna["y.linked"] <- list(NULL)
  r <- infer_quiet(kna, system = "xy", seed = 100)
  expect_true(is.na(r$agreed.sex[1])) # was '*F'
  expect_true(is.na(r$x.linked.sex[1]))
  expect_equal(r$agreed.sex[2], "M")
})

test_that("duplicated individual names keep their assignments (change 7)", {
  kd <- k_xy
  nm <- indNames(kd$x.linked)
  nm[2] <- nm[1]
  for (s in c("y.linked", "x.linked", "gametolog")) indNames(kd[[s]]) <- nm
  r <- infer_quiet(kd, system = "xy", seed = 100)
  ref <- infer_quiet(k_xy, system = "xy", seed = 100)
  expect_identical(r$id, nm)
  expect_identical(r$x.linked.sex, ref$x.linked.sex) # was all NA
})

test_that("verbosity (change 6)", {
  expect_silent(infer_quiet(k_xy, system = "xy", seed = 1))
  out <- capture.output(gl.infer.sex(k_xy, system = "xy", seed = 1,
                                     verbose = 2))
  expect_true(any(grepl("Completed: gl.infer.sex", out)))
  expect_true(any(grepl("gametologs \\(at least 5", out)))
})

test_that("k-means always splits the sample into two sexes", {
  fem <- which(LBP@other$ind.metrics$sex == "F")
  kf <- k_xy
  kf$x.linked <- k_xy$x.linked[fem, ]
  kf["y.linked"] <- list(NULL)
  kf["gametolog"] <- list(NULL)
  r <- infer_quiet(kf, system = "xy", seed = 100)
  expect_equal(as.vector(table(r$agreed.sex)), c(94L, 68L))
})

test_that("the caller's random number stream is restored (change 5)", {
  set.seed(1)
  a <- runif(1)
  set.seed(1)
  invisible(infer_quiet(k_xy, system = "xy", seed = 100))
  expect_equal(runif(1), a) # was different
  set.seed(1)
  invisible(infer_quiet(k_xy, system = "xy"))
  expect_equal(runif(1), a)
})

test_that("platypus, xy", {
  r <- infer_quiet(keep_quiet(platypus.gl, "xy"), system = "xy", seed = 1)
  expect_equal(as.vector(table(r$agreed.sex)), c(54L, 27L))
})
