#' Advanced Options for Mixed-effect modeling functions.
#'
#' Change default settings for gremlin models.
#'
#' @aliases gremlinControl
#' @param cctol Convergence criteria tolerances (Meyer 2007, 2019).
#' @param ezero Effective zero to be used, values less than this number are
#'   treated as zero and fixed to this value.
#' @param einf Effective infinite value to be used, values are limited to a
#'   to this variable as a maximum.
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
#'   Meyer, K. 2007. WOMBAT - a tool for mixed model analyses in quantitative
#'   genetics by restricted maximum likelihood (REML). Journal of Zhejiang
#'   University SCIENCE B 8(11):815-821.
#'
#'   Meyer, K. 2019. WOMBAT A program for mixed model analyses by restricted
#'   maximum likelihood. User Notes. 27 September 2019.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#'   str(gremlinControl())
#'
#' @export
gremlinControl <- function(cctol = c(5e-4, 1e-8, 1e-3, NULL),
	ezero = 1e-8, einf = 1e30, step = 0.3, lambda = TRUE,
	algorithm = NULL, algArgs = list()){

  stopifnot(is.list(algArgs))

 return(structure(list(cctol = cctol, ezero = ezero, einf = einf, step = step,
		lambda = lambda,
		algorithm = algorithm, algArgs = algArgs),
	class = c("gremlinControl")))
}  #<-- end `gremlinControl()`

