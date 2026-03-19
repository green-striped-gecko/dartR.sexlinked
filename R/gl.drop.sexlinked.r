#'@name gl.drop.sexlinked
#'@title Used to remove loci that are sex linked
#'@description
#' This function is deprecated. Please use \code{\link{gl.filter.sexlinked}} 
#' instead.
#'
#' @param x Name of the genlight object containing the SNP data. This genlight
#' object needs to contain the sex of the individuals [required].
#' @param system String that declares the sex-determination system of the 
#' species: 'zw' or 'xy' [required].
#' 
#' @details
#' This function has been deprecated and replaced by gl.filter.sexlinked in 
#' order to keep consistency with other functions (gl.report -> gl.filter).
#'
#' @return A warning.
#' 
#' @author Custodian: Diana Robledo-Ruiz -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @references
#' \itemize{
#' \item Robledo‐Ruiz, D. A., Austin, L., Amos, J. N., Castrejón‐Figueroa, J.,
#'  Harley, D. K., Magrath, M. J., Sunnucks, P., & Pavlova, A. (2023). 
#'  Easy‐to‐use R functions to separate reduced‐representation genomic datasets
#'   into sex‐linked and autosomal loci, and conduct sex assignment. Molecular 
#'   Ecology Resources, 00, 1-21.
#'  }
#' 
#' @export
#'
gl.drop.sexlinked <- function(x, 
                              system = NULL) {
  warning("`gl.drop.sexlinked()` is deprecated; please use `gl.filter.sexlinked()`
          instead.", 
          call. = FALSE)
}