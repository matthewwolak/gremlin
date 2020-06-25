################  Extract `logLik` and `AIC`/`AICc`    #############
# Create S3methods that use generic in package `stats`
######   LOG-LIKELIHOOD   ######
#' Methods to extract log-likelihood and information criterion of a gremlin
#' model.
#'
#' Extracts the log-likelihood or AIC from a gremlin model fit.
#' 
#' Function \code{npar.gremlin} returns an object with attributes \code{n.fxd}
#' and \code{n.bndry} which give additional information about the parameters
#' estimated and contributing to the overall \code{df} of the model. \code{n.fxd}
#' returns the total number of parameters (No. fixed effects + No. (co)variance
#' components) minus the number of parameters constrained to a certain value. Thus,
#' \code{n.fxd} represents the number of parameters that can vary and, as a 
#' consequence, affect the log-likelihood.
#'
#' The attribute \code{n.bndry} reports the number of parameters that were
#' restrained to stay inside the boundaries of allowable parameter space (e.g.,
#' a variance that was not allowed to be negative).
#'
#' @aliases logLik.gremlin npar.gremlin AIC.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param \dots Additional arguments.
#' @param k A numeric value for the penalty per parameter. Default is 2, as in
#'   classic AIC.
#' @param fxdDf A logical indicating whether to penalize according to the number
#'   of fixed effect parameters. Since only models fit by REML can be compared,
#'   these must always be the same and so become a constant. Hence, the default
#'   is \code{FALSE}.
#'
#' @return \code{numeric} values for the log-likelihood, the number of
#'   parameters estimated by the model (sum of fixed effects and random effect
#'   (co)variance components), and Akaike's Information Criterion.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' grS <- gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#' logLik(grS)
#' AIC(grS)
#' @export
#' @importFrom stats logLik
logLik.gremlin <- function(object, ...){
  val <- object$itMat[nrow(object$itMat), "loglik"]
  structure(val,
	    nobs = nobs.gremlin(object),
#	    nall = nobs.gremlin(object),  #<-- as in `lme4::logLik.merMod()`
  	    df = npar.gremlin(object),
            class = "logLik")
}


######  Number of parameters estimated ######
#' @rdname logLik.gremlin
#' @export
npar.gremlin <- function(object){
  n <- object$grMod$modMats$nb + object$grMod$p
  con <- object$grMod$conv
  fxd <- sum(con == "F")  #<-- number of parameters fixed in the model
  bndry <- sum(con == "B")  #<-- number of parameters restrained to a boundary

  structure(n,                    # number of parameters of the model
	    n.fxd = n - fxd,      # number of 'free' parameters
            n.bndry = bndry)      # DIFFERENT: number of parameters at a boundary
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
  if(fxdDf) p <- df else p <- df - object$grMod$modMats$nb  
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
#' grS <- gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#' nobs(grS)
#' @export
#' @importFrom stats nobs
#TODO check that calculates correct value if NAs are dropped: see `modMats()`
nobs.gremlin <- function(object, use.fallback = FALSE, ...){
  object$grMod$modMats$ny
}








############          Time to run model    ################
#' Time to execute the gremlin model
#'
#' Extract the length of time to fit the model.
#' 
#' @aliases runtime
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param \dots Further arguments to be passed to the methods.
#'
#' @return A numeric of class \sQuote{difftime} with an attribute of units
#' (e.g., seconds or minutes).
#' @author \email{matthewwolak@@gmail.com}
runtime <- function(object, ...){
  attr(object, "endTime") - attr(object, "startTime")
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
#' mod11 <- gremlin(WWG11 ~ sex - 1,
#'	random = ~ sire,
#'	data = Mrode11)
#' mod11red <- gremlinR(WWG11 ~ sex - 1, data = Mrode11)
#' anova(mod11, mod11red)
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

  grMods <- as.logical(vapply(dots, FUN = is, FUN.VALUE = NA, "gremlin"))
  if(sum(c(is(object, "gremlin"), grMods)) < 2){
    stop("At least 2 models must be of class 'gremlin' (`is(object, 'gremlin')`)")
  }
  if(any(!c(is(object, "gremlin"), grMods))){
    stop("All models must be of class 'gremlin' (`is(object, 'gremlin')`)")
  }

  if(is(object, "gremlin")){
    mods <- c(list(object), dots[grMods])
  } else mods <- dots[grMods]

  nobs.vec <- vapply(mods, nobs, 1L)
  if(var(nobs.vec) > 0){
    stop("Models were not all fitted to the same size of dataset")
  }
  if(is.null(mNms <- model.names)){
    mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, grMods)], FUN = .safeDeparse, "")  
  }

  if(any(duplicated(mNms)) || max(nchar(mNms)) > 200){
    warning("Failed to find model names, assigning generic names")
    mNms <- paste0("model", seq_along(mNms))
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
  calls <- lapply(mods, FUN = function(x) getCall(x))
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
    "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower.tail = FALSE), #TODO introduce mixtures
    row.names = names(mods), check.names = FALSE)
  class(tabout) <- c("anova", class(tabout))
  forms <- lapply(lapply(calls, FUN = "[[", "random"), FUN = deparse)
  rcovForms <- lapply(lapply(calls, FUN = "[[", "rcov"), FUN = deparse)
  structure(tabout,
    heading = c(header, "Models:",
      paste(paste(rep.int(names(mods), lengths(forms)), unlist(forms), sep = ": "),
        unlist(rcovForms), sep = " + [ rcov ] ")))
 
}









################################################################################
################  Extract Fixed Effect Estimates    #############
# Create S3methods that use generic in package `stats`
# adapted from lme4::fixef()
#' Fixed Effect Estimates of \code{class} \sQuote{gremlin}
#'
#' Extracts the fixed effect estimates from a model of \code{class} \sQuote{gremlin}.
#' 
#' @aliases fixef fixed.effects fixef.gremlin
#' @param object An object of \code{class} \sQuote{gremlin}.
#' @param add.dropped A \code{logical} value indicating whether fixed effects dropped by
#'   gremlin, due to rank deficiencies in the fixed effect design matrix, should
#'   be included with \code{NA} values.
#' @param \dots Additional arguments.
#'
#' @return A \code{numeric} vector of fixed effect estimates.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' fixef(gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11))
#' @importFrom nlme fixef
#' @export fixef
#' @method fixef gremlin
#' @export
fixef.gremlin <- function(object, add.dropped = FALSE, ...){

  if(add.dropped){
    warning(cat("add.dropped = TRUE not implemented\n"))
  }

 structure(object$grMod$sln[1:object$grMod$modMats$nb],
           names = dimnames(object$grMod$modMats$X)[[2L]])
}


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
#' grS <- gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#' residuals(grS)
#' @export
#' @importFrom stats naresid residuals
# type = c("working", "response", "deviance", "pearson", "partial")
residuals.gremlin <- function(object,
	type = "response",
	scaled = FALSE, ...){
  if(missing(type)) type <- "response" else type <- match.arg(type)
  r <- as.vector(object$grMod$r)
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
#TODO how to make `scaled` if multivariate (>1 residual variance/covariance)?
  if(scaled) res <- res / sqrt(object$grMod$thetav[object$grMod$p])
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
#'     \item{formulae }{Function call and model fixed, random, and residual
#'       formulae.}
#'     \item{runtime }{A \code{numeric} of class \sQuote{difftime} containing
#'       the length of time to run the model. See how this is handled in
#'       \code{\link{update.gremlin}}.}
#'     \item{lambda }{A \code{logical} indicating if the model was transformed
#'       to the variance ratio, or \code{lambda} scale.}
#'     \item{residQuants }{A named \code{vector} listing summary output for the
#'       model residuals.}
#'     \item{varcompSummary }{Table of variance components and approximate
#'       standard errors (calculated from the inverse of the average information
#'       matrix). If a (co)variance component is fixed or at the boundary of
#'       its parameter space then an \code{NA} is returned for the standard error
#'       and a column with constraint types is added to the table. Alternative
#'       methods (e.g., profile likelihood CIs) should be pursued for obtaining
#'       uncertainties associated with fixed or boundary parameters.}
#'     \item{varcompSampCor }{A \code{matrix} containing the sampling correlations
#'       of the (co)variance components. Note this is on the underlying \code{nu}
#'       scale that the model is fitting.}
#'     \item{coefficients }{Table of fixed effects and standard errors (calculated
#'       from the corresponding diagonal elements of the inverse of the
#'       coefficient matrix, transformed where necessary).}
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{gremlin}}
#' @examples
#' grS <- gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#' summary(grS)
#' @export
#' @importFrom stats quantile
#' @method summary gremlin
summary.gremlin <- function(object, ...){
  nit <- nrow(object$itMat)
  nvc <- object$grMod$p      # No. of (co)variance parameters
  formulae <- list(fn = object$grMod$call[[1L]],
	fxd = object$grMod$call[["formula"]],
	random = object$grMod$call[["random"]],
        rcov = object$grMod$call[["rcov"]])

  ########
  residQuants <- quantile(residuals(object, "response", scaled = TRUE), na.rm=TRUE)

  ########
  vcItMatInd <- grep("_theta", dimnames(object$itMat)[[2L]])
  vcItMatNames <- sapply(strsplit(dimnames(object$itMat)[[2L]][vcItMatInd],
    split = "_"), FUN = "[", i = 1)
  varcompSummary <- cbind("Estimate" = object$itMat[nit, vcItMatInd, drop = TRUE],
		"Std. Error" = NA)
    dimnames(varcompSummary)[[1L]] <- vcItMatNames
    if(!is.null(object$grMod$AI)){
      invAI <- try(solve(object$grMod$AI), silent = TRUE)
      if(inherits(invAI, "try-error")){
        varcompSummary[, "Std. Error"] <- rep(NA, nrow(varcompSummary))
          attr(varcompSummary, "condition") <- attr(invAI, "condition")$message
        varcompSampCor <- matrix(NA, nrow = nvc, ncol = nvc)
          dimnames(varcompSampCor) <- list(dimnames(varcompSummary)[[1L]],
                                           dimnames(varcompSummary)[[1L]])
      } else{
          if(object$grMod$lambda){
            varcompSummary[, "Std. Error"] <- sqrt(nuVar2thetaVar_lambda(object))
          } else varcompSummary[, "Std. Error"] <- sqrt(diag(invAI))
          varcompSampCor <- cov2cor(invAI)
        }
    } else{
        varcompSampCor <- matrix(NA, nrow = nvc, ncol = nvc)
          dimnames(varcompSampCor) <- list(dimnames(varcompSummary)[[1L]],
                                           dimnames(varcompSummary)[[1L]])
      }
    # Remove SEs of fixed or boundary parameters
    if(any(FBparams <- object$grMod$conv %in% c("F", "B"))){
      varcompSummary[which(FBparams), "Std. Error"] <- NA
      conCol <- rep("", nrow(varcompSummary))
        conCol[which(FBparams)] <- sapply(as.character(object$grMod$conv[FBparams]),
              FUN = switch, "F" = "Fixed", "B" = "Boundary")
      varcompSummary <- as.data.frame(varcompSummary)
      varcompSummary$Constraint <- conCol
     }

  ########
  coefVarScale <- ifelse(object$grMod$lambda, object$grMod$sigma2e, 1.0)
  fxdSummary <- cbind(as(fixef(object), "matrix"),
        sqrt(coefVarScale * object$grMod$Cinv_ii[1:object$grMod$modMats$nb]))
    #TODO consider reporting `|z value|` instead
    fxdSummary <- cbind(fxdSummary, fxdSummary[, 1] / fxdSummary[, 2])
    colnames(fxdSummary) <- c("Estimate", "Std. Error", "z value")
	dimnames(fxdSummary)[[1L]] <- names(fixef(object))

 return(structure(list(logLik = logLik(object),
		formulae = formulae,
		runtime = runtime(object),
		lambda = object$grMod$lambda,
		residQuants = residQuants,
		varcompSummary = varcompSummary,
		varcompSampCor = varcompSampCor,
		coefficients = fxdSummary),
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
  cat("\n Linear mixed model fit by REML ['", x$formula$fn, "']")
  cat("\n REML log-likelihood:", round(x$logLik, digits+2), "\n")
    cat("\t lambda:", x$lambda, "\n")
  cat("\n elapsed time for model:", round(x$runtime, digits), "\n")

  # Adapted from `lme4::.prt.resids`
    cat("\n Scaled residuals:\n")
    resids <- setNames(zapsmall(x$residQuants, digits + 1L),
	c("Min", "1Q", "Median", "3Q", "Max"))
    print(resids, digits = digits, ...)

  cat("\n (co)variance parameters:", paste(as.expression(x$formulae$random)), "\n")
  cat("\t\t    rcov:", paste(as.expression(x$formulae$rcov)), "\n")
    print(as.data.frame(x$varcompSummary), digits = digits, ...)

  cat("\n (co)variance parameter sampling correlations:\n")
    print(as.data.frame(x$varcompSampCor), digits = digits, ...)
  if(!is.null(conOut <- attr(x$varcompSummary, "condition"))){
    cat("(co)variance parameter Std. Errors/sampling correlations not available\n",
      "Failed to invert AI matrix\n", conOut, "\n")
  }

  cat("\n Fixed effects:", paste(as.expression(x$formulae$fxd)), "\n")
    #See `printCoefmat()`
    print(as.data.frame(x$coefficients), digits = digits, ...)

}









