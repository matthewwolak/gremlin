################  Extract `logLik` and `AIC`/`AICc`    #############
# Create S3methods that use generic in package `stats`
######   LOG-LIKELIHOOD   ######
#' Methods to extract log-likelihood and information criterion of a gremlin
#' model.
#'
#' Extracts the log-likelihood or AIC from a gremlin model fit.
#' 
#' @aliases logLik.gremlin AIC.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param \dots Additional arguments.
#' @param k A numeric value for the penalty per parameter. Default is 2, as in
#'   classic AIC.
#' @param fxdDf A logical indicating whether to penalize according to the number
#'   of fixed effect parameters. Since only models fit by REML can be compared,
#'   these must always be the same and so become a constant. Hence, the default
#'   is \code{FALSE}.
#'
#' @return A \code{numeric} value for either the log-likelihood and the number of
#'   parameters estimated by the model (sum of fixed effects and random effect
#'   (co)variance components) or Akaike's Information Criterion.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' mod11 <- gremlinR(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	Gstart = matrix(0.1), Rstart = matrix(0.4),
#'	maxit = 50, v = 2, algit = "EM")
#' logLik(mod11)
#' AIC(mod11)
#' @export
#' @importFrom stats logLik
logLik.gremlin <- function(object, ...){
  val <- object$itMat[nrow(object$itMat), "loglik"]
  #TODO attr(val, "nall") <- object$
  attr(val, "df") <- object$modMats$nb + nrow(object$theta)
  class(val) <- "logLik"
 val
}


######   AIC   ######
#' @rdname logLik.gremlin
#' @export
#' @importFrom stats AIC
AIC.gremlin <- function(object, ..., k = 2, fxdDf = FALSE){
  ll <- logLik(object)
    df <- attr(ll, "df")
  #TODO determine which is better: asreml (not including fixed df) or lme4 (includes)
  ## might depend on REML vs ML
  if(fxdDf) p <- df else p <- df - object$modMats$nb  
  #TODO implement conditional AIC (df=n-1 levels of each random effect)
  ## I think below does what `lme4::lmer()` implements
  ###TODO Will need to accommodate covariance parameters (count lower triangles)
  val <- -2*ll[[1L]] + k*p
  #TODO document type = c("mAIC") attribute and reference (see GLMM faq for refs)
  attr(val, "type") <- "mAIC"
  class(val) <- "AIC"
 val
}
#TODO BIC
## See how asreml calculates this (equation 2.15): -2ll + t * log(residual df)
#TODO make an extractor for all 3 (like `lme4::llikAIC()`)











############          Number of observations in fitted model    ################
#' Number of observations in data from gremlin model fit objects
#'
#' Extract the number of 'observations' in a gremlin model fit.
#' 
#' @aliases nobs.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param use.fallback logical: should fallback methods be used to try to guess
#'   the value? Included for compatibility.
#' @param \dots Further arguments to be passed to the methods.
#'
#' @return A single number, usually an \code{integer}, but can be \code{NA}.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' mod11 <- gremlinR(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	Gstart = matrix(0.1), Rstart = matrix(0.4),
#'	maxit = 50, v = 2, algit = "EM")
#' nobs(mod11)
#' @export
#' @importFrom stats nobs
#TODO check that calculates correct value if NAs are dropped: see `modMats()`
nobs.gremlin <- function(object, use.fallback = FALSE, ...){
  object$modMats$ny
}








############          REML Likelihood Ratio Tests using `anova()`    ##########
#' anova() for gremlin objects
#'
#' REML Likelihood Ratio Tests for gremlin models using anova()
#' 
#' @aliases anova.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param \dots Additional objects of \code{class} \sQuote{gremlin}.
#' @param model.names Optional character vector with model names to be used in
#'   the anova table
#'
#' @return A \code{data.frame} containing the nested comparison of model
#'   \code{object}s via a REML likelihood ratio test.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' mod11 <- gremlinR(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	Gstart = matrix(0.1), Rstart = matrix(0.4),
#'	maxit = 50, v = 2, algit = "EM")
#' logLik(mod11)
#' @export
#' @importFrom stats anova getCall pchisq
#adapted from `lme4::anovaLmer()`
anova.gremlin <- function(object, ..., model.names = NULL){
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  # Define temporary/internal function to unlist an `lapply()`
  .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
  # Define temporary/internal functions: taken from `lme4`
    .safeDeparse <- function(x, collapse = " "){
      paste(deparse(x, 500L), collapse = collapse)
    }
    .abbrDeparse <- function(x, width = 60){
      r <- deparse(x, width)
      if(length(r) > 1) paste(r[1], "...") else r
    }
  #TODO reconcile likelihoods with `lm()` where REML = TRUE so can compare
  ##TODO add `|` below to allow lm
  grMods <- as.logical(vapply(dots, FUN = is, FUN.VALUE = NA, "gremlin"))
  if(sum(c(is(object) == "gremlin", grMods)) < 2){
    stop("At least 2 models must be of class 'gremlin' (`is(object) == gremlin`)")
  }

  if(is(object) == "gremlin"){
    mods <- c(list(object), dots[grMods])
  } else mods <- dots[grMods]

  nobs.vec <- vapply(mods, nobs, 1L)
  if(var(nobs.vec) > 0)
    stop("Models were not all fitted to the same size of dataset")

  if(is.null(mNms <- model.names)){
    mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, grMods)], FUN = .safeDeparse, "")  
  }

  if(any(duplicated(mNms)) || max(nchar(mNms)) > 200){
    warning("Failed to find model names, assigning generic names")
    nNms <- paste0("model", seq_along(mNms))
  }

  if(length(mNms) != length(mods))
    stop("Number of models and 'model.names' vector are not of the same length")

  names(mods) <- sub("@env$", '', mNms)
  ###########
  lls <- lapply(mods, logLik)
  ## Order models by increasing degrees of freedom
  llo <- order(Df <- vapply(lls, FUN = attr, FUN.VALUE = numeric(1), "df"))
  mods <- mods[llo]
  lls <- lls[llo]
  Df <- Df[llo]
  calls <- lapply(mods, getCall)
  data <- lapply(calls, FUN = "[[", "data")
  if(!all(vapply(data, FUN = identical, FUN.VALUE = NA, data[[1]]))){
    stop("All models must be fit to the same data object")
  }

  fxdForms <- lapply(calls, FUN = "[[", "formula")
  if(!all(vapply(fxdForms, FUN = identical, FUN.VALUE = NA, fxdForms[[1]]))){
    badFxd <- which(!vapply(fxdForms, FUN = identical, FUN.VALUE = NA, fxdForms[[1]]))
    cat("Different fixed effect formulas for model(s)", names(mods)[badFxd], ":\n")
    stop(cat("All models must use the same fixed effects formula (as compared to model '", names(mods)[1], "')\n"))
  }

  header <- paste("Data:", .abbrDeparse(data[[1]]))
  header <- c(header, paste("Fixed formula:", deparse(fxdForms[[1]])))
  subset <- lapply(calls, FUN = "[[", "subset")
  if(!all(vapply(subset, FUN = identical, FUN.VALUE = NA, subset[[1]]))){
    stop("All models must use the same subset")
  }
  if(!is.null(subset[[1]])){
    header <- c(header, paste("Subset:", .abbrDeparse(subset[[1]])))
  }
  
  llk <- unlist(lls)
  chisq <- 2 * pmax(0, c(NA, diff(llk)))
  dfChisq <- c(NA, diff(Df))
  tabout <- data.frame(Df = Df,
    AIC = .sapply(lls, AIC),
#TODO    BIC = .sapply(lls, BIC),
    logLik = llk,
    deviance = -2*llk,
    Chisq = chisq,
    "Chi Df" = dfChisq,
    "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower.tail = FALSE),
    row.names = names(mods), check.names = FALSE)
  class(tabout) <- c("anova", class(tabout))
  forms <- lapply(lapply(calls, FUN = "[[", "random"), FUN = deparse)
  structure(tabout,
    heading = c(header, "Models:",
      paste(rep.int(names(mods), lengths(forms)), unlist(forms), sep = ": ")))
 
}









################################################################################
################  Extract Residuals    #############
# Create S3methods that use generic in package `stats`
#' Residuals of \code{class} \sQuote{gremlin}
#'
#' Residuals of \code{class} \sQuote{gremlin}.
#' 
#' @aliases residuals.gremlin resid.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param type The type of residuals which should be returned. Only implement
#'   \dQuote{response} currently. Can be abbreviated.
#' @param scaled Logical value indicating whether to scale residuals by the
#'   residual standard deviation.
#' @param \dots Additional arguments.
#'
#' @return A \code{numeric} vector of residuals.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' mod11 <- gremlinR(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	Gstart = matrix(0.1), Rstart = matrix(0.4),
#'	maxit = 50, v = 2, algit = "EM")
#' residuals(mod11)
#' @export
#' @importFrom stats naresid residuals
# type = c("working", "response", "deviance", "pearson", "partial")
residuals.gremlin <- function(object,
	type = "response",
	scaled = FALSE, ...){
  if(missing(type)) type <- "response" else type <- match.arg(type)
  r <- object$residuals
#  res <- switch(type, 
#	working = ,
#	response = r,
#	deviance = ,
#	pearson = ,  #<-- shouldn't allow weights - that is what R-structure does
#	partial = r)
  res <- r  #<-- since don't use `switch()`
  res <- naresid(object$na.action, res)
#TODO make a `predict.gremlin()` method
  if(type == "partial" | type == "pearson"){
    stop("partial and pearson residuals not implemented")
#  res <- res + predict(object, type = "terms")
  }
  if(scaled) res <- res / sqrt(object$theta[object$modMats$nG+1])
 res
}
#TODO######   hatvalues (see lme4) & predict.gremlin    ############





################################################################################
#' Gremlin model summary.
#'
#' Summarize and print results of linear mixed model fitted with gremlin.
#'
#' @aliases summary.gremlin print.summary.gremlin
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
#' @export
#' @importFrom stats quantile
#' @method summary gremlin
summary.gremlin <- function(object, ...){
  nit <- nrow(object$itMat)
  nvc <- nrow(object$dLdtheta)     # No. of (co)variance components
  formulae <- list(fxd = object$call[["formula"]],
	random = object$call[["random"]]) #FIXME need to include R (or combine with G)

  ########
  residQuants <- quantile(residuals(object, "response", scaled = TRUE), na.rm=TRUE)

  ########
  varcompSummary <- cbind("Estimate" = object$itMat[nit, 1:nvc, drop = TRUE],
		"Std. Error" = NA)
    dimnames(varcompSummary)[[1L]] <- dimnames(object$itMat[nit, 1:nvc, drop = FALSE])[[2L]]
    if(!is.null(object$AI)){
      invAI <- solve(object$AI)
      varcompSummary[, "Std. Error"] <- sqrt(diag(invAI))
      varcompSampCor <- cov2cor(invAI)
    } else{
        varcompSampCor <- matrix(NA, nrow = nvc, ncol = nvc)
          dimnames(varcompSampCor) <- list(dimnames(varcompSummary)[[1L]], dimnames(varcompSummary)[[1L]])
      }

  ########
  fxdSummary <- object$sln[1:object$modMats$nb, , drop = FALSE]
    fxdSummary[, 2] <- sqrt(fxdSummary[, 2])
    #TODO consider reporting `|z value|` instead
    fxdSummary <- cbind(fxdSummary, fxdSummary[, 1] / fxdSummary[, 2])
    colnames(fxdSummary) <- c("Solution", "Std. Error", "z value")
	dimnames(fxdSummary)[[1L]] <- object$modMats$X@Dimnames[[2L]]

 return(structure(list(logLik = logLik(object),
		formulae = formulae,
		residQuants = residQuants,
		varcompSummary = varcompSummary,
		varcompSampCor = varcompSampCor,
		fxdSummary = fxdSummary),
	class = c("summary.gremlin", "list")))
} 











################################################################################
#' @method print summary.gremlin
#' @rdname summary.gremlin
#' @export
#' @importFrom stats setNames
print.summary.gremlin <- function(x, 
	digits = max(3, getOption("digits") - 3), ...){
#TODO calculate convergence criteria & print if REML converged
## also print if parameters changed by >XX%
  cat("\nLinear mixed model fit by REML ['gremlin']")
  cat("\nREML log-likelihood:", round(x$logLik, digits), "\n")

  # Adapted from `lme4::.prt.resids`
    cat("\nScaled residuals:\n")
    resids <- setNames(zapsmall(x$residQuants, digits + 1L),
	c("Min", "1Q", "Median", "3Q", "Max"))
    print(resids, digits = digits, ...)

  cat("\nRandom effects:", paste(as.expression(x$formulae$random)), "\n")
    print(as.data.frame(x$varcompSummary), digits = digits, ...)

  cat("\nRandom effect sampling correlations:\n")
    print(as.data.frame(x$varcompSampCor), digits = digits, ...)

  cat("\nFixed effects:", paste(as.expression(x$formulae$fxd)), "\n")
    #See `printCoefmat()`
    print(as.data.frame(x$fxdSummary), digits = digits, ...)

}









