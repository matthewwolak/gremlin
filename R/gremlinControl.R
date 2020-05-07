#' Advanced Options for Mixed-effect modeling functions.
#'
#' Change default settings for gremlin models.
#'
#' @aliases gremlinControl
#' @param cctol Convergence criteria tolerances. TODO
#' @param ezero Effective zero to be used, values less than this number are
#'   treated as zero and fixed to this value.
#' @param step A \code{numeric} value for scaling the proposed parameter updates.
#' @param lambda A \code{logical} indicating whether a residual variance should
#'   be factored out of the mixed model equations.
#' @param algorithm A \code{character} naming the function to use to decide
#'   subsequent parameters in the REML iterations.
#' @param algArgs A \code{list} of function arguments to be given to functions
#'   named in the \code{algorithm} argument.
#'
#' @return A \code{list} of class \code{gremlinControl} to be used by
#'   \code{gremlinSetup} and later functions when fitting the model.
#' @references
#'   Meyer. Convergence Checks
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#'   str(gremlinControl())
#'
#' @export
gremlinControl <- function(cctol = c(5e-4, 1e-8, 1e-3, NULL),
	ezero = 1e-8, step = 0.3, lambda = TRUE,
	algorithm = NULL, algArgs = list()){

  stopifnot(is.list(algArgs))

 return(structure(list(cctol = cctol, ezero = ezero, step = step,
		lambda = lambda,
		algorithm = algorithm, algArgs = algArgs),
	class = c("gremlinControl")))
}  #<-- end `gremlinControl()`

