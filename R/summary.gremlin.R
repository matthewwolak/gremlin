################  Extract `logLik` and `AIC`/`AICc`    #############
# Create S3methods that use generic in package `stats`
######   LOG-LIKELIHOOD   ######
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
print.summary.gremlin <- function(x, digits = max(3, getOption("digits") - 3), ...){
#TODO calculate convergence criteria and print whether REML converged & if parameters changed by >XX%
  cat("\n log-likelihood:", round(x$logLik, digits))
  cat("\n Variance components:", paste(as.expression(x$formulae$random)), "\n\n")
  print(as.data.frame(x$varcompSummary), digits = digits, ...)

  cat("\n Fixed effects:", paste(as.expression(x$formulae$fxd)), "\n\n")
  #See `printCoefmat()`
  print(as.data.frame(x$fxdSummary), digits = digits, ...)

}




