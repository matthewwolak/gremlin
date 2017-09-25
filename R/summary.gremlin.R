################  Extract `logLik` and `AIC`/`AICc`    #############
# Create S3methods that use generic in package `stats`
######   LOG-LIKELIHOOD   ######
#' Methods to extract log-likelihood and information criterion of a gremlin
#' model.
#'
#' Extracts the log-likelihood from a gremlin model fit.
#' 
#' @aliases logLik.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param \dots Additional arguments to be passed to control the model fitting.
#'
#' @return A \code{numeric} value for the log-likelihood.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' mod11 <- gremlinR(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	Gstart = matrix(0.1), Rstart = matrix(0.4),
#'	maxit = 50, v = 2, algit = "EM")
#' logLik(mod11)
#' @export
#' @import stats
logLik.gremlin <- function(object, ...){
  val <- object$itMat[nrow(object$itMat), "loglik"]
  #TODO attr(val, "nall") <- object$
  #TODO attr(val, "df") <- object$
  class(val) <- "logLik"
 val
}
#TODO######   AIC    ############
#TODO?######   BIC    ############






################################################################################
#' Gremlin model summary.
#'
#' Summarize and print results of linear mixed model fitted with gremlin.
#'
#' @aliases summary.gremlin print.summary.gremlin
#' @export
#' @param object,x An object of \code{class} \sQuote{gremlin} or
#'   \sQuote{summary.gremlin}.
#' @param digits An \code{integer} used for number formatting with
#'   \sQuote{signif()}. 
#' @param \dots Additional arguments to be passed to control the output.
#'
#' @return A \code{list} of class \code{summary.gremlin} or a printed value
#'   to the screen with no return values.
#'   \describe{
#'     \item{logLik }{Model log-likelihood.}
#'     \item{formulae }{Model fixed, random, and residual formulae.}
#'     \item{varcompSummary }{Table of variance components and approximate
#'       standard errors (calculated from the inverse of the average information
#'       matrix).}
#'     \item{fxdSummary }{Table of fixed effects and standard errors (calculated
#'       from the corresponding diagonal elements of the inverse of the
#'       coefficient matrix).}
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{gremlin}}
summary.gremlin <- function(object, ...){
  nit <- nrow(object$itMat)

  formulae <- list(fxd = NULL, random = NULL) #FIXME need to combine G and R

  varcompSummary <- cbind(Est = object$itMat[nit, 1:nrow(object$dLdtheta)],
		SE = sqrt(diag(solve(object$AI))))
  fxdSummary <- object$sln[1:object$modMats$nb, ]
    fxdSummary[, 2] <- sqrt(fxdSummary[, 2])
    colnames(fxdSummary)[2L] <- "SE"
	dimnames(fxdSummary)[[1L]] <- object$modMats$X@Dimnames[[2L]]

 return(structure(list(logLik = logLik(object),
		formulae = formulae,
		varcompSummary = varcompSummary,
		fxdSummary = fxdSummary),
	class = c("summary.gremlin", "list")))
} 












################################################################################
#' @method print summary.gremlin
#' @rdname summary.gremlin
#' @export
print.summary.gremlin <- function(x, 
	digits = max(3, getOption("digits") - 3), ...){
#TODO calculate convergence criteria & print if REML converged
## also print if parameters changed by >XX%
  cat("\n log-likelihood:", round(x$logLik, digits))
  cat("\n Variance components:", paste(as.expression(x$formulae$random)), "\n\n")
  print(as.data.frame(x$varcompSummary), digits = digits, ...)

  cat("\n Fixed effects:", paste(as.expression(x$formulae$fxd)), "\n\n")
  #See `printCoefmat()`
  print(as.data.frame(x$fxdSummary), digits = digits, ...)

}




