#'@name gl.drop.sexlinked
#'@title Removes loci that are sex linked (deprecated)
#'@description
#' This function is deprecated. Please use \code{\link{gl.filter.sexlinked}}
#' instead. It warns and then runs \code{gl.filter.sexlinked} with the same
#' arguments, so scripts written for it keep working.
#'
#' @param x Name of the genlight object containing the SNP data. This genlight
#' object needs to contain the sex of the individuals [required].
#' @param system String that declares the sex-determination system of the
#' species: 'zw' or 'xy' [required].
#' @param ... Further arguments passed to \code{gl.filter.sexlinked}
#' (ncores, plot.display, plot.theme, plot.colors, plot.file, plot.dir,
#' verbose).
#'
#' @details
#' This function has been deprecated and replaced by gl.filter.sexlinked in
#' order to keep consistency with other functions (gl.report -> gl.filter).
#'
#' @return What \code{gl.filter.sexlinked} returns: a genlight object with the
#' autosomal loci, or NULL.
#'
#' @author Author(s): Diana Robledo-Ruiz. Custodian: Diana Robledo-Ruiz --
#'   Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @references
#' \itemize{
#' \item Robledo‐Ruiz, D. A., Austin, L., Amos, J. N., Castrejon‐Figueroa, J.,
#'  Harley, D. K., Magrath, M. J., Sunnucks, P., & Pavlova, A. (2023). 
#'  Easy‐to‐use R functions to separate reduced‐representation genomic datasets
#'   into sex‐linked and autosomal loci, and conduct sex assignment. Molecular 
#'   Ecology Resources, 00, 1-21.
#'  }
#'
#' @seealso \code{\link{gl.filter.sexlinked}}
#'
#' @export
#'
gl.drop.sexlinked <- function(x,
                              system = NULL,
                              ...) {
  .Deprecated("gl.filter.sexlinked", package = "dartR.sexlinked")
  # Re-evaluate the caller's own call under the new name, so the history
  # entry gl.filter.sexlinked() records names the caller's object and
  # arguments (x = LBP, system = "xy"), as a direct call would
  cl <- match.call()
  cl[[1]] <- quote(gl.filter.sexlinked)
  eval(cl, parent.frame())
}
