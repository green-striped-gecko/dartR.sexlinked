#' @name gl.keep.sexlinked
#' @title Keeps loci that are sex linked
#' @family data manipulation
#' @description
#' This function identifies sex-linked and autosomal loci present in a SNP
#' dataset (genlight object) using individuals with known sex. It identifies
#' five types of loci: w-linked or y-linked, sex-biased, z-linked or
#' x-linked, gametologous and autosomal.
#'
#' This function returns a list with 5 elements, including one dataframe and
#' 4 genlight objects with sex-linked loci, and displays 4 plots.
#'
#' @param x Name of the genlight object containing the SNP data. This genlight
#' object needs to contain the sex of the individuals. See explanation in
#' details [required].
#' @param system String that declares the sex-determination system of the
#' species: 'zw' or 'xy' [required].
#' @param ncores Number of processes to be used in parallel operation. If ncores
#' > 1 parallel operation is activated [default 1].
#' @param plot.display If TRUE, displays the four output plots. See
#' explanation in details [default TRUE].
#' @param plot.theme Theme for the plots [default theme_dartR()].
#' @param plot.colors Not implemented; the plots use fixed colours for each
#' category of loci [default NULL].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' The genlight object must contain in \code{gl@other$ind.metrics} a column
#' named "sex" (case-insensitive) in which individuals with known sex are
#' assigned 'F' or 'Female' for females, and 'M' or 'Male' for males
#' (case-insensitive). The rows of \code{ind.metrics} must be in the same
#' order as the individuals in the genlight object. The function ignores
#' individuals that are assigned anything else or nothing at all
#' (unknown sex). At least one female and one male are required. SilicoDArT
#' (presence/absence) data are not accepted.
#'
#' Loci are classified as in \code{\link{gl.report.sexlinked}}, which
#' describes the tests and thresholds used.
#'
#' The display of plots can be turned off (\code{plot.display = FALSE}), but
#' we strongly encourage you to always inspect the output plots at least once
#' to make sure everything is working properly.
#'
#'\strong{ Function's output }
#'
#' This function returns a list of 5 elements: \itemize{
#' \item {$results.table > Table with statistics (columns) for each loci
#' (rows), as returned by \code{\link{gl.report.sexlinked}}}
#' \item {$w.linked or $y.linked > Genlight object with w-linked/y-linked loci}
#' \item {$sex.biased > Genlight object with sex-biased scoring rate loci}
#' \item {$z.linked or $x.linked > Genlight object with z-linked/x-linked loci}
#' \item {$gametolog > Genlight object with gametologs}
#' }
#' A category with no loci is returned as NULL. Each genlight object keeps all
#' individuals, carries the locus metrics of its own loci only, and records
#' this call in its history.
#'
#' And displays four plots:\itemize{
#' \item {A BEFORE plot based on loci call rate by sex, with w/y-linked loci 
#' colored in yellow and sex-biased loci in blue}
#' \item {An AFTER plot based on loci call rate by sex, with only sex-linked
#'  loci}
#' \item {A BEFORE plot based on loci heterozygosity by sex, with z/x-linked 
#' loci colored in orange and gametologs in green}
#' \item {An AFTER plot based on loci heterozygosity by sex, with only 
#' sex-linked loci}
#' }
#' The plots are not returned; use \code{plot.file} to save them.
#'
#' @return A list of 5 elements (see Details).
#'
#' @author Author(s): Diana Robledo-Ruiz. Custodian: Diana Robledo-Ruiz --
#'   Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' LBP_sexLinked <- gl.keep.sexlinked(x = LBP, system = "xy", 
#' plot.display = TRUE, ncores = 1)
#' LBP_sexLinked$gametolog
#'
#' @references
#' \itemize{
#' \item Robledo-Ruiz, D. A., Austin, L., Amos, J. N., Castrejon-Figueroa, J.,
#'  Harley, D. K., Magrath, M. J., Sunnucks, P., & Pavlova, A. (2023).
#'  Easy to use R functions to separate reduced representation genomic datasets
#'   into sex linked and autosomal loci, and conduct sex assignment. Molecular
#'   Ecology Resources, 00, 1-21.
#'  }
#'
#' @importFrom stats chisq.test
#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
#' @importFrom foreach foreach "%dopar%"
#'
#' @export

gl.keep.sexlinked <- function(x,
                              system = NULL,
                              ncores = 1,
                              plot.display = TRUE,
                              plot.theme = theme_dartR(),
                              plot.colors = NULL,
                              plot.file = NULL,
                              plot.dir = NULL,
                              verbose = NULL) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.display <- FALSE
  }
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # SET COLOURS #not yet implemented...
  if (is.null(plot.colors)) {
    plot.colors <- c("#2171B5", "#6BAED6")
  } else {
    if (length(plot.colors) > 2) {
      if (verbose >= 2) {
        cat(warn(
          "  More than 2 colors specified, only the first 2 are used\n"
          ))
      }
      plot.colors <- plot.colors[1:2]
    }
  }
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  # Heterozygosity is undefined for presence/absence data, so SilicoDArT is
  # not accepted
  datatype <- utils.check.datatype(x,
                                 accept = c("dartR",
                                              "genlight", 
                                              "SNP"),
                                   verbose = verbose)
  
  if (is.null(system)) {
    stop(
      error(
        "You must specify the sex-determination system with the parameter 'system' ('zw' or 'xy')."
      )
    )
  } else {
    if (!(system == 'zw' | system == 'xy')) {
      stop(error(
        "Parameter 'system' must be 'zw' or 'xy'."
        ))
    }
  }
  
  # Extract sex per individual
  metrics <- x@other$ind.metrics
  
  # Locate the "sex" column regardless of case
  sex.cols <- grep("^sex$",
                   names(metrics),
                   ignore.case = TRUE,
                   value = TRUE)
  
  if (length(sex.cols) == 0) {
    stop(error("Could not find any column named 'sex' (case-insensitive)."))
  }
  
  if (length(sex.cols) > 1) {
    cat(warn(sprintf(
      "  Multiple columns matched 'sex' (case-insensitive): %s. Using the first: '%s'.",
      paste(sex.cols, collapse = ", "),
      sex.cols[1]
    )))
  }
  
  sex.col <- sex.cols[1]
  
  # ind.metrics rows track individuals 1:1, so sex is read by row position
  # rather than by matching an 'id' column against indNames(x)
  if (nrow(metrics) != nInd(x)) {
    stop(error(paste0(
      "The number of rows in ind.metrics (", nrow(metrics),
      ") does not match the number of individuals (", nInd(x), ").\n"
    )))
  }
  
  # Pull values and force upper case
  sex.values <- toupper(metrics[[sex.col]])
  is.F <- sex.values %in% c("F", "FEMALE")
  is.M <- sex.values %in% c("M", "MALE")
  
  if (verbose > 1)
    message(report(paste(
      "  Detected ",
      sum(is.F),
      " females and ",
      sum(is.M),
      " males.",
      sep = ""
    )))
  
  # Both sexes are needed to compare call rate and heterozygosity
  if (sum(is.F) == 0 || sum(is.M) == 0) {
    stop(error(paste0(
      "Found ", sum(is.F), " females and ", sum(is.M), " males; at least ",
      "one of each is needed. Females and males in the 'sex' column must be ",
      "'F' and 'M', or 'FEMALE' and 'MALE' (case-insensitive).\n"
    )))
  }
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    # Stop the workers even if the function exits with an error
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
  }
  
  # Transform genotypes to matrix (loci in rows) and subset by sex;
  # drop = FALSE keeps a matrix when a sex has a single individual
  gen <- t(as.matrix(x))
  gen.F <- gen[, is.F, drop = FALSE]
  gen.M <- gen[, is.M, drop = FALSE]
  
  # Names of the columns that depend on the sex-determination system
  if (system == "zw") {
    col.hemi <- "w.linked"
    col.dip  <- "z.linked"
  } else {
    col.hemi <- "y.linked"
    col.dip  <- "x.linked"
  }
  
  # Test for independence of sex and a 2x2 table of counts, given in the
  # order F row first column, M row first column, F row second column, M row
  # second column. Returns the estimate and its p-value.
  sex.test <- function(counts) {
    obs <- matrix(counts, nrow = 2, ncol = 2)
    
    # See if it's possible to use chisq test
    if (sum(obs) >= 1000) {
      # Convert zeros to 1 so that chisq.test does not return NaN when a
      # row or column sums to zero
      obs[obs == 0] <- 1
      res <- chisq.test(obs, correct = FALSE)
      c(unname(res$statistic), res$p.value)
    } else {
      # Fisher's exact test (because there are observations with less
      # than 5). It accepts zeros, so the observed counts are used as is
      res <- fisher.test(obs)
      c(unname(res$estimate), res$p.value)
    }
  }
  
  # Apply sex.test to every row of a matrix of counts
  run.tests <- function(counts) {
    # Bind the foreach loop variable to avoid an R CMD check note
    i <- NULL
    if (ncores > 1) {
      res <- foreach::foreach(i = seq_len(nrow(counts)),
                              .combine = rbind,
                              .export = "sex.test") %dopar% {
        sex.test(counts[i, ])
      }
      matrix(res, ncol = 2)
    } else {
      t(vapply(seq_len(nrow(counts)),
               function(i) sex.test(counts[i, ]),
               numeric(2)))
    }
  }
  
  ##################### 1. Sex-linked loci by scoring rate
  
  # Create a results table with loci as row names and index number per locus
  table <- data.frame(index = c(1:nrow(gen)), row.names = row.names(gen))
  
  # Count missing (NA) and add as column to table
  table$count.F.miss <- rowSums(is.na(gen.F))
  table$count.M.miss <- rowSums(is.na(gen.M))
  
  # Count scored ("0", "1" or "2") and add as column to table
  table$count.F.scored <- rowSums(!is.na(gen.F))
  table$count.M.scored <- rowSums(!is.na(gen.M))
  
  if (verbose > 1) {
    if (ncores > 1) {
      message(report("  Starting phase 1. Working in parallel..."))
    } else {
      message(report("  Starting phase 1. May take a while..."))
    }
  }
  
  # Test for independence of sex and missingness
  res <- run.tests(as.matrix(table[, c("count.F.miss", "count.M.miss",
                                       "count.F.scored", "count.M.scored")]))
  table$ratio   <- res[, 1]
  table$p.value <- res[, 2]
  
  scoringRate.F <- scoringRate.M <- heterozygosity.F <- heterozygosity.M <- NA
  # Adjust p-values for multiple comparisons (False discovery rate)
  table$p.adjusted <- p.adjust(table$p.value, method = "fdr")
  
  # Calculate scoring rate for females and males and add to results table
  table$scoringRate.F <- table$count.F.scored / (table$count.F.scored +
                                                   table$count.F.miss)
  
  table$scoringRate.M <- table$count.M.scored / (table$count.M.scored +
                                                   table$count.M.miss)
  
  ##### 1.1 W-linked or Y-linked loci
  # Loci scored in at most 10% of the homogametic sex (males for zw,
  # females for xy) with a significant sex effect on call rate
  if (system == "zw") {
    rate.absent <- table$scoringRate.M
  } else {
    rate.absent <- table$scoringRate.F
  }
  sig.miss <- table$p.adjusted <= 0.01
  table[[col.hemi]] <- rate.absent <= 0.1 & sig.miss
  
  ##### 1.2 Loci with sex-biased scoring rate
  table$sex.biased <- sig.miss & !table[[col.hemi]]
  
  table.hemi      <- table[table[[col.hemi]], ]
  table.sexbiased <- table[table$sex.biased, ]
  
  ##### 1.3 Plot BEFORE vs AFTER
  if (verbose > 1) {
    message(report(
      "  Building call rate plots."
      ))
  }
  
  table.autosomal <- table[!table[[col.hemi]] & !table$sex.biased, ]
  
  BEF.mis <- ggplot2::ggplot(table.autosomal, 
                             aes(x = scoringRate.F, y = scoringRate.M)) +
    geom_point(color = 'grey33') +
    geom_point(data = table.sexbiased, color = 'dodgerblue3') +
    geom_point(data = table.hemi, color = 'gold') +
    ggtitle("BEFORE keeping only sex-linked loci") +
    xlab("Call rate Females") +
    ylab("Call rate Males") +
    xlim(0, 1) + ylim(0, 1) +
    plot.theme
  
  AFT.mis <- ggplot2::ggplot(table.sexbiased, 
                             aes(x = scoringRate.F, y = scoringRate.M)) +
    geom_point(color = 'dodgerblue3') +
    geom_point(data = table.hemi, color = 'gold') +
    ggtitle("AFTER keeping only sex-linked loci") +
    xlab("Call rate Females") +
    ylab("Call rate Males") +
    xlim(0, 1) + ylim(0, 1) +
    plot.theme
  
  #################### 2. Sex-linked loci by heterozygosity
  # Count heterozygotes ("1") and add as column to results table
  table$count.F.het <- rowSums(gen.F == 1, na.rm = TRUE)
  table$count.M.het <- rowSums(gen.M == 1, na.rm = TRUE)
  
  # Count homozygotes ("0" or "2") and add as column to results table
  table$count.F.hom <- rowSums(
    gen.F != 1,
    # Ignores NAs
    na.rm = TRUE
  )
  table$count.M.hom <- rowSums(
    gen.M != 1,
    # Ignores NAs
    na.rm = TRUE
  )
  
  if (verbose > 1){
    message(report(
      "  Starting phase 2. May take a while..."
      ))
  }
  
  # Apply test for independence of sex and heterozygosity, excluding
  # w/y-linked loci and loci with sex-biased score
  table$stat         <- NA_real_
  table$stat.p.value <- NA_real_
  tested <- which(!table[[col.hemi]] & !table$sex.biased)
  if (length(tested) > 0) {
    res <- run.tests(as.matrix(table[tested, c("count.F.het", "count.M.het",
                                               "count.F.hom", "count.M.hom")]))
    table$stat[tested]         <- res[, 1]
    table$stat.p.value[tested] <- res[, 2]
  }
  
  # Adjust p-values for multiple comparisons (False discovery rate, 
  # least conservative)
  table$stat.p.adjusted <- p.adjust(table$stat.p.value, method = "fdr")
  
  # Calculate for heterozygosity per sex and add to results table
  table$heterozygosity.F <- table$count.F.het / 
    (table$count.F.het + table$count.F.hom)
  
  table$heterozygosity.M <- table$count.M.het / 
    (table$count.M.het + table$count.M.hom)
  
  ##### 2.1 Z-linked or X-linked loci AND gametologs
  # Among loci with a significant sex effect on heterozygosity, those more
  # heterozygous in the homogametic sex (males for zw, females for xy) are
  # z/x-linked; the rest are gametologs
  sig.het <- !is.na(table$stat.p.adjusted) & table$stat.p.adjusted <= 0.01
  if (system == "zw") {
    higher <- table$heterozygosity.M > table$heterozygosity.F
  } else {
    higher <- table$heterozygosity.F > table$heterozygosity.M
  }
  table[[col.dip]]  <- sig.het & higher %in% TRUE
  table$gametolog   <- sig.het & higher %in% FALSE
  
  table.dip     <- table[table[[col.dip]], ]
  table.gametol <- table[table$gametolog, ]
  
  ##### 2.2 Plot BEFORE vs AFTER
  if (verbose > 1) {
    message(report(
      "  Building heterozygosity plots."
      ))
  }
  
  is.autosomal <- !table[[col.hemi]] & !table$sex.biased &
    !table[[col.dip]] & !table$gametolog
  table.autosomal <- table[is.autosomal, ]
  
  BEF.het <- ggplot2::ggplot(table.autosomal,
                             aes(x = heterozygosity.F, y = heterozygosity.M)) +
    geom_point(color = 'grey33') +
    geom_point(data = table.gametol, color = 'chartreuse3') +
    geom_point(data = table.dip, color = 'darkorange1') +
    ggtitle("BEFORE keeping only sex-linked loci") +
    xlab("% Heterozygous Females") +
    ylab("% Heterozygous Males") +
    xlim(0, 1) + ylim(0, 1) +
    plot.theme
  
  AFT.het <- ggplot2::ggplot(table.gametol,
                             aes(x = heterozygosity.F, y = heterozygosity.M)) +
    geom_point(color = 'chartreuse3') +
    geom_point(data = table.dip, color = 'darkorange1') +
    ggtitle("AFTER keeping only sex-linked loci") +
    xlab("% Heterozygous Females") +
    ylab("% Heterozygous Males") +
    xlim(0, 1) + ylim(0, 1) +
    plot.theme
  
  if (verbose > 1) {
    message(report(
      "  Done building heterozygosity plots."
      ))
  }
  
  #################### 3. Create output of function
  ##### 3.1 Save the indices of each category of loci to later subset x
  a <- table[table[[col.hemi]], "index"]
  b <- table[table$sex.biased, "index"]
  c <- table[table[[col.dip]], "index"]
  d <- table[table$gametolog, "index"]
  
  if (system == "zw") {
    lab.hemi <- "W-linked"
    lab.dip  <- "Z-linked"
  } else {
    lab.hemi <- "Y-linked"
    lab.dip  <- "X-linked"
  }
  
  if (verbose > 1) message("**FINISHED** \nTotal of analyzed loci: ", nrow(table), ".\n",
          "Kept ", length(a) + length(b) + length(c) + length(d), " sex-linked loci:\n",
          "   ",    length(a), " ", lab.hemi, " loci (yellow)\n",
          "   ",    length(b), " sex-biased loci (blue)\n",
          "   ",    length(c), " ", lab.dip, " loci (orange)\n",
          "   ",    length(d), " gametologs (green).\n",
          "And removed ", sum(is.autosomal), " autosomal loci (grey).")
  
  ##### 3.2 Subset x object
  # Returns NULL for an empty category. loc.metrics are re-subset from x
  # because the genlight '[' method (unlike dartR's) leaves them untouched
  this.call <- match.call()
  keep.loci <- function(idx) {
    if (length(idx) == 0) {
      return(NULL)
    }
    y <- x[, idx]  # Loci are columns
    if (!is.null(x@other$loc.metrics)) {
      y@other$loc.metrics <- x@other$loc.metrics[idx, , drop = FALSE]
    }
    # ADD TO HISTORY
    nh <- length(y@other$history)
    y@other$history[[nh + 1]] <- this.call
    y
  }
  
  A <- keep.loci(a)
  B <- keep.loci(b)
  C <- keep.loci(c)
  D <- keep.loci(d)
  
  #################### 4. Output
  rlist <- list(table, A, B, C, D)
  names(rlist) <- c("results.table", col.hemi, "sex.biased", col.dip,
                    "gametolog")
  
  p4 <- BEF.mis + AFT.mis + BEF.het + AFT.het
  
  if (plot.display) {
    print(p4)
  }
  # Optionally save the plot ---------------------
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p4,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # FLAG SCRIPT END ---------------
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  # ----------------------
  
  # RETURN
  return(rlist)
}
