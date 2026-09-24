#' @name gl.infer.sex
#' @title Uses sex-linked loci to infer sex of individuals
#' @family unmatched report
#' @description
#' This function uses the output of function gl.keep.sexlinked (list of 5
#' elements) to infer the sex of all individuals in the dataset.
#' It uses 3 types of sex-linked loci (W-/Y-linked, Z-/X-linked, and
#' gametologs), assigns a preliminary genetic sex for each type of sex-linked
#' loci available, and outputs an agreed sex.
#'
#' This function produces as output a dataframe with individuals in rows and 11
#' columns.
#'
#' @param gl_sexlinked The output of function gl.keep.sexlinked (complete
#' list with 5 elements). See explanation in "Details" section [required].
#' @param system String that declares the sex-determination system of the
#' species: 'zw' or 'xy'. It must be the same system used in
#' gl.keep.sexlinked [required].
#' @param seed User-defined integer for repeatability purposes. If not provided
#' by user, it is chosen randomly. See "Details" section [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' Parameter \code{gl_sexlinked} must be the name of the output object (a
#' list of 5 elements) produced by function \code{gl.keep.sexlinked}, run with
#' the same \code{system}. Parameter \code{seed} must be an integer that will
#' be used on the KMeans algorithm used by the function. If \code{seed} is not
#' provided, the value chosen is not reported, so we highly recommend choosing
#' the seed to guarantee repeatability. The random number stream of the R
#' session is restored when the function ends.
#'
#' Each type of sex-linked loci is used as follows:
#' \itemize{
#' \item {W-linked or Y-linked loci: each individual is assigned to the
#' heterogametic sex if it has more called than missing genotypes, and to the
#' homogametic sex otherwise.}
#' \item {Z-linked or X-linked loci: used only if there are at least 2 loci.
#' Individuals are split into two groups by KMeans on their proportions of
#' heterozygous and homozygous genotypes.}
#' \item {Gametologs: used only if there are at least 5 loci. Only the 5
#' gametologs with the smallest adjusted p-value in the results table are
#' used, and they are split into two groups by KMeans as for Z-/X-linked
#' loci.}
#' }
#' KMeans always returns two groups, so the Z-/X-linked and gametolog
#' assignments assume that both sexes are present in the dataset. If only one
#' sex is present, part of the individuals will still be assigned to the
#' other sex.
#'
#' Note that this function was created with the explicit intent that a human
#' checks the evidence for the sex assignments that do NOT agree for all
#' types of sex-linked loci (called "indefinite sex assignments" and denoted
#' as "*M", "*F" or "*?" in the last column of dataframe output). This human
#' can then use their criterion to validate these assignments.
#'
#'\strong{ Function's output }
#'
#' This function creates a dataframe with one row per individual and 11
#' columns:
#' \itemize{
#' \item {id > Individuals' ID.}
#' \item {w.linked.sex or y.linked.sex > Sex inferred using w-linked or y-linked
#' loci.}
#' \item {#called > Number of W-linked or Y-linked loci for which the individual
#' had a called genotype (cf. missing genotype).}
#' \item {#missing > Number of W-linked or Y-linked loci for which the
#' individual had a missing genotype (cf. called genotype).}
#' \item {z.linked.sex or x.linked.sex > Sex inferred using z-linked or x-linked
#' loci.}
#' \item {#Het.z or #Het.x > Number of z-linked or x-linked loci for which the
#' individual is heterozygous.}
#' \item {#Hom.z or #Hom.x > Number of z-linked or x-linked loci for which the
#' individual is homozygous.}
#' \item {gametolog.sex > Sex inferred using gametologs.}
#' \item {#Het.g > Number of the 5 gametologs used for which the individual is
#' heterozygous.}
#' \item {#Hom.g > Number of the 5 gametologs used for which the individual is
#' homozygous.}
#' \item {agreed.sex > Agreed sex: 'F' or 'M' if all available preliminary
#' sex-assignments match (i.e., definite sex assignment); '*F' or '*M' if NOT
#' all preliminary sex-assignments match (i.e., indefinite sex assignment),
#' following the majority; '*?' if there is no majority; and NA if no
#' preliminary sex-assignment is available.}
#' }
#' Columns for a type of sex-linked loci that could not be used are NA.
#'
#' @return A dataframe with one row per individual and the 11 columns
#' described in "Details".
#' @author Author(s): Diana Robledo-Ruiz. Custodian: Diana Robledo-Ruiz --
#'   Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' LBP_sexLinked <- gl.keep.sexlinked(x = LBP, system = "xy", 
#' plot.display = FALSE, ncores = 1)
#' inferred.sexes <- gl.infer.sex(gl_sexlinked = LBP_sexLinked, system = "xy", 
#' seed = 100)
#' inferred.sexes
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
#' @seealso \code{\link{gl.keep.sexlinked}}
#'
#' @importFrom stats kmeans
#' @importFrom stats na.omit
#'
#' @export

gl.infer.sex <- function(gl_sexlinked,
                         system = NULL,
                         seed = NULL,
                         verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # Parameters check
  if (is.null(system)) {
    stop(error(
      "You must specify the sex-determination system with the parameter 'system' ('zw' or 'xy')."
    ))
  } else {
    if (!(system == 'zw' | system == 'xy')) {
      stop(error(
        "Parameter 'system' must be 'zw' or 'xy'."
        ))
    }
  }
  
  if (system == "xy") {
    name1 <- "y.linked"
    name2 <- "x.linked"
  } else {
    name1 <- "w.linked"
    name2 <- "z.linked"
  }
  
  if (!is.list(gl_sexlinked) ||
      is.null(gl_sexlinked[["results.table"]])) {
    stop(error(
      "Parameter 'gl_sexlinked' must be the list returned by gl.keep.sexlinked().\n"
    ))
  }
  
  # A subset set to NULL is still a named element, so a missing name means
  # that the list was built with the other system
  if (!all(c(name1, name2) %in% names(gl_sexlinked))) {
    other <- if (system == "xy") "zw" else "xy"
    stop(error(
      paste0("'gl_sexlinked' has no '", name1, "' and '", name2,
             "' elements: it was built with system = '", other,
             "'. Use the same system in gl.keep.sexlinked() and ",
             "gl.infer.sex().\n")
    ))
  }
  
  gl1 <- gl_sexlinked[[name1]]
  gl2 <- gl_sexlinked[[name2]]
  
  # Gametologs
  gl3    <- gl_sexlinked[["gametolog"]]
  table  <- gl_sexlinked$results.table  # Retrieve table
  all    <- table[table$gametolog == TRUE, ]
  all    <- all[order(all$stat.p.adjusted), ]  # Order from smallest p-value
  useful <- row.names(all[1:5, ])              # Keep name of only top 5 gametologs
  
  obj <- Filter(Negate(is.null), list(gl1, gl2, gl3))
  if (length(obj) == 0) {
    stop(error(
      "No sex-linked loci (W-/Y-linked, Z-/X-linked or gametologs) found in 'gl_sexlinked'.\n"
    ))
  }
  ind_names <- indNames(obj[[1]])
  n_ind     <- length(ind_names)
  if (any(sapply(obj, nInd) != n_ind)) {
    stop(error(
      "The genlight objects in 'gl_sexlinked' do not have the same individuals.\n"
    ))
  }
  
  # kmeans() uses the random number stream; restore the caller's stream on
  # exit so that setting 'seed' here does not change it
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
            add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }
  
  # Random seed if not specified by user
  if (is.null(seed)) {
    seed <- sample.int(65535, 1)
  }
  
  # Make sex assignment per type of sex-linked loci (Functions declared below)
  # W/Y-linked
  if (!is.null(gl1)) {
    w <- W.sex(gl1, system = system)
    if (verbose >= 2) {
      cat(report("  Assigning sexes based on W-linked/Y-linked loci.\n"))
    }
  } else {
    if (verbose >= 2) {
      cat(warn(
        "  No W-linked/Y-linked loci were found. Assigning NAs to sex",
        "assignment based on W-linked/Y-linked loci.\n"
      ))
    }
    w <- data.frame(
      W.sex = rep(NA, n_ind),
      n0.w  = rep(NA, n_ind),
      n1.w  = rep(NA, n_ind)
    )
  }
  # Z/X-linked
  if (!is.null(gl2) && gl2@n.loc >= 2) {
    z <- Z.sex(gl2, system = system, seed = seed)
    if (verbose >= 2) {
      cat(report("  Assigning sexes based on Z-linked/X-linked loci.\n"))
    }
  } else {
    if (verbose >= 2) {
      cat(warn(
        "  Not enough Z-linked/X-linked loci (at least 2 are needed).",
        "Assigning NAs to sex assignment based on Z-linked/X-linked loci.\n"
      ))
    }
    z <- data.frame(
      Z.sex = rep(NA, n_ind),
      n1.z  = rep(NA, n_ind),
      n0.z  = rep(NA, n_ind)
    )
  }
  # Gametologs
  if (!is.null(gl3) && gl3@n.loc >= 5) {
    g <- g.sex(gl3,
               system = system,
               seed = seed,
               useful = useful)
    if (verbose >= 2) {
      cat(report("  Assigning sexes based on gametologous loci.\n"))
    }
  } else {
    if (verbose >= 2) {
      cat(warn(
        "  Not enough gametologs (at least 5 are needed). Assigning NAs to",
        "sex assignment based on gametologous loci.\n"
      ))
    }
    g <- data.frame(
      g.sex = rep(NA, n_ind),
      n1.g  = rep(NA, n_ind),
      n0.g  = rep(NA, n_ind)
    )
  }
  
  # Put them all together (rows are in the order of the individuals)
  A <- data.frame(w, z, g)
  
  # Function to conciliate assignments
  Fun <- function(x, y, z) {
    d <- as.vector(na.omit(c(x, y, z)))
    if (length(d) == 0) {
      return(NA_character_)  # No assignment available
    }
    if (length(unique(d)) == 1) {
      return(d[1])
    }
    n_f <- sum(d == "F")
    n_m <- sum(d == "M")
    if (n_f == n_m) {
      return("*?")  # No majority
    }
    if (n_f > n_m) "*F" else "*M"
  }
  
  # Add last column of agreed sexes
  A$agreed.sex <- mapply(Fun, A$W.sex, A$Z.sex, A$g.sex, USE.NAMES = FALSE)
  
  if (system == 'xy') {
    names <- c(
      'y.linked.sex',
      '#called',
      '#missing',
      'x.linked.sex',
      '#Het.x',
      '#Hom.x',
      'gametolog.sex',
      '#Het.g',
      '#Hom.g',
      'agreed.sex'
    )
  } else {
    names <- c(
      'w.linked.sex',
      '#called',
      '#missing',
      'z.linked.sex',
      '#Het.z',
      '#Hom.z',
      'gametolog.sex',
      '#Het.g',
      '#Hom.g',
      'agreed.sex'
    )
  }
  
  colnames(A) <- names
  
  A <- cbind(id = ind_names, A)
  
  # Row names must be unique, so they are kept only when the IDs are
  if (!anyDuplicated(ind_names)) {
    row.names(A) <- ind_names
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(A)
}

############################### 1. W.sex function
### Map NAs (missing) and scored (called) to 1Dim in [-1,1], IF x<0, F else M

W.sex <- function(gl, system = NULL) {
  w <- as.matrix(gl)
  w[is.na(w)] <- 3
  
  n0.w <- unname(rowSums(w == 0 | w == 2 | w == 1, na.rm = TRUE))
  n1.w <- unname(rowSums(w == 3, na.rm = TRUE))
  
  # Calculate proportion
  sex.score <- function(f, m) {
    return((-f + m) / (f + m))
  }
  
  c2 <- sex.score(n0.w, n1.w)
  
  if (system == 'xy') {
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }
  
  W.sex <- ifelse(c2 < 0, lab0, lab1)
  
  Y <- data.frame(W.sex, n0.w, n1.w)
  return(Y)
}

############################### 2. Z.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Hom

Z.sex <- function(gl, system = NULL, seed = 42) {
  z <- as.matrix(gl)
  
  n0.z <- unname(rowSums(z == 0 | z == 2, na.rm = TRUE))
  n1.z <- unname(rowSums(z == 1, na.rm = TRUE))
  
  Z  <- cbind(n0.z, n1.z) / (n0.z + n1.z)
  ok <- !is.na(Z[, 1])  # Individuals with no called genotype are left out
  
  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z[ok, , drop = FALSE], 2)
  
  # Keep clusters in the order of the individuals (NA if left out)
  cluster     <- rep(NA_integer_, length(n0.z))
  cluster[ok] <- km$cluster
  
  label <- cluster[which.max(n1.z)] # Largest number of '1'
  
  if (system == 'xy') {
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }
  
  # Assign sex
  Z.sex <- ifelse(cluster == label, lab1, lab0)
  
  Y <- data.frame(Z.sex, n1.z, n0.z)
  return(Y)
}

############################### 3. ZWg.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Het

g.sex <-  function(gl,
                   system = NULL,
                   seed = 42,
                   useful = useful) {
  z <- as.matrix(gl)
  z <- z[, useful]
  
  n0.g <- unname(rowSums(z == 0 | z == 2, na.rm = TRUE))
  n1.g <- unname(rowSums(z == 1, na.rm = TRUE))
  
  Z  <- cbind(n0.g, n1.g) / (n0.g + n1.g)
  ok <- !is.na(Z[, 1])  # Individuals with no called genotype are left out
  
  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z[ok, , drop = FALSE], 2)
  
  # Keep clusters in the order of the individuals (NA if left out)
  cluster     <- rep(NA_integer_, length(n0.g))
  cluster[ok] <- km$cluster
  
  label <- cluster[which.max(n1.g)] # Largest number of '1'
  
  if (system == 'xy') {
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }
  
  # Assign sex (HERE IS THE OPPOSITE)
  g.sex <- ifelse(cluster == label, lab0, lab1)
  
  Y <- data.frame(g.sex, n1.g, n0.g)
  return(Y)
}
