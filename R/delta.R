################################################################################
# Delta Method methods
#######################

##############
# Generic
##############
#' Delta Method to Calculate Standard Errors for Functions of (Co)variances.
#'
#' Calculates the standard error for results of simple mathematical functions of
#'   (co)variance parameters using the delta method.
#'
#' The delta method (e.g., Lynch and Walsh 1998, Appendix 1; Ver Hoef 2012) uses
#'   a Taylor series expansion to approximate the moments of a function of
#'   parameters. Here, a second-order Taylor series expansion is implemented to
#'   approximate the standard error for a function of (co)variance parameters.
#'   Partial first derivatives of the function are calculated by algorithmic
#'   differentiation with \code{\link[stats]{deriv}}.
#'
#' Though \code{deltaSE} can calculate standard errors for non-linear functions
#'   of (co)variance parameters from a fitted \code{gremlin} model, it is limited
#'   to non-linear functions constructed by mathematical operations such as the
#'   arithmetic operators \code{+}, \code{-}, \code{*}, \code{/} and \code{^},
#'   and single-variable functions such as  \code{exp} and \code{log}. See 
#'   \code{\link[stats]{deriv}} for more information.
#'
#' @aliases deltaSE deltaSE.default deltaSE.formula deltaSE.list
#' @param calc A character \code{expression}, \code{formula}, or list (of
#'   \code{formula} or \code{expression}) expressing the mathematical function
#'   of (co)variance component for which to calculate standard errors.
#' @param object A fitted model object of \code{class} \sQuote{gremlin}.
#' @param scale A \code{character} indicating whether to calculate the function
#'   and standard error on the original data scale (\dQuote{theta}) or
#'   on the underlying scale to which (co)variance components are transformed
#'   for the model fitting calculations (\dQuote{nu}). Defaults to
#'   \dQuote{theta} if not specified.
#' @return A \code{data.frame} containing the \dQuote{Estimate} and
#'   \dQuote{Std. Error} for the mathematical function(s) of (co)variance
#'   components.
#' @references 
#'   Lynch, M. and B. Walsh 1998. Genetics and Analysis of Quantitative Traits.
#'   Sinauer Associates, Inc., Sunderland, MA, USA.
#'
#'   Ver Hoef, J.M. 2012. Who invented the delta method? The American
#'   Statistician 66:124-127. DOI: 10.1080/00031305.2012.687494
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[stats]{deriv}}
#' @importFrom stats deriv
#' @examples
#'   # Calculate the sum of the variance components 
#'     grS <- gremlin(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#'     deltaSE(Vsum ~ V1 + V2, grS)
#'     deltaSE("V1 + V2", grS)  #<-- alternative
#'
#'   # Calculate standard deviations (with standard errors) from variances
#'     ## Uses a `list` as the first (`calc`) argument
#'     ### All 3 below: different formats to calculate the same values
#'     deltaSE(list(SD1 ~ sqrt(V1), SDresid ~ sqrt(V2)), grS)  #<-- formulas
#'     deltaSE(list(SD1 ~ sqrt(G.sire), SDresid ~ sqrt(ResVar1)), grS) 
#'     deltaSE(list("sqrt(V1)", "sqrt(V2)"), grS)  #<-- list of character expressions
#'
#'   # Additive Genetic Variance calculated from observed Sire Variance
#'     ## First simulate Full-sib data
#'     set.seed(359)
#'     noff <- 5     #<-- number of offspring in each full-sib family
#'     ns <- 100     #<-- number of sires/full-sib families
#'     VA <- 1       #<-- additive genetic variance
#'     VR <- 1       #<-- residual variance
#'     datFS <- data.frame(id = paste0("o", seq(ns*noff)),
#'       sire = rep(paste0("s", seq(ns)), each = noff))
#'     ## simulate mid-parent breeding value (i.e., average of sire and dam BV)
#'     ### mid-parent breeding value = 0.5 BV_sire + 0.5 BV_dam
#'     #### var(mid-parent BV) = 0.25 var(BV_sire) + 0.25 var(BV_dam) = 0.5 var(BV) 
#'     datFS$midParBV <- rep(rnorm(ns, 0, sqrt(0.5*VA)), each = noff)
#'     ## add to this a Mendelian sampling deviation to get each offspring BV
#'     datFS$BV <- rnorm(nrow(datFS), datFS$midParBV, sqrt(0.5*VA))
#'     datFS$r <- rnorm(nrow(datFS), 0, sqrt(VR))  #<-- residual deviation
#'     datFS$pheno <- rowSums(datFS[, c("BV", "r")]) 
#'     # Analyze with a sire model
#'     grFS <- gremlin(pheno ~ 1, random = ~ sire, data = datFS)
#'     # calculate VA as 2 times the full-sib/sire variance
#'     deltaSE(VAest ~ 2*V1, grFS)
#'     # compare to expected value and simulated value
#'     VA  #<-- expected
#'     var(datFS$BV)  #<-- simulated (includes Monte Carlo error)
#'
#'   # Example with `deltaSE(..., scale = "nu")
#'   ## use to demonstrate alternative way to do same calculation of inverse
#'   ## Average Information matrix of theta scale parameters when lambda = TRUE
#'   ### what is done inside gremlin::nuVar2thetaVar_lambda 
#'     dOut <- deltaSE(thetaV1 ~ V1*V2, grS, "nu")  #<-- V2 is sigma2e
#'     aiFnOut <- nuVar2thetaVar_lambda(grS)[1]  #<-- variance (do sqrt below)
#'     stopifnot(abs(dOut[, "Std. Error"] - sqrt(aiFnOut)) < 1e-10)
#'
#' @export
deltaSE <- function(calc, object, scale = c("theta", "nu")){
  UseMethod("deltaSE", calc)
}



##############
# Default
##############
#' @describeIn deltaSE Default method
#' @export
#' @importFrom stats as.formula
deltaSE.default <- function(calc, object, scale = c("theta", "nu")){

  if(!inherits(object, "gremlin")){
    stop(cat("Must supply an object of class", dQuote(gremlin), "\n"))
  }

  cl <- match.call()
  cl_calc <- cl[[match("calc", names(cl))]]
  if(inherits(cl_calc, "character")){
    fmla <- as.formula(paste("~", eval(calc)))
  } else{
      warning("calc must be a character expression")
      fmla <- as.formula(paste("~", eval(cl_calc)))
    }

 deltaSE.formula(fmla, object, scale)
}  #<-- end deltaSE.default





##############
# Formula
##############
#' @describeIn deltaSE Formula method
#' @export
#' @importFrom stats deriv
deltaSE.formula <- function(calc, object, scale = c("theta", "nu")){

  if(!inherits(object, "gremlin")){
    stop(cat("Must supply an object of class", dQuote(gremlin), "\n"))
  }

  RHS <- calc[[length(calc)]]

  # determine if theta or nu scale calculations to be performed
  if(missing(scale)) sc <- "theta" else sc <- match.arg(scale)
  if(sc == "theta"){
    if(object$grMod$lambda){
      invAI <- nuAI2thetaAIinv_lambda(object)
    } else invAI <- solve(object$grMod$AI)  #<-- TODO if non-lambda nu!=theta
    covlist <- as.list(object$grMod$thetav)
   
  } else{  #<-- if "nu"
      invAI <- solve(object$grMod$AI)
      covlist <- as.list(matlist2vech(object$grMod$nu))
      if(object$grMod$lambda) covlist[which(covlist == 1.0)] <- object$grMod$sigma2e
    }  #<-- end if/else theta/nu 


  # check whether names of covlist (theta/nu) are used in `RHS`
  ## first create generic "V" names
  covlistVnms <- paste0("V", seq(length(covlist)))
  onames <- any(names(covlist) %in% all.vars(RHS))
  vnames <- any(covlistVnms %in% all.vars(RHS))
  if(!onames && !vnames){
    stop(cat("Variables in calc are not names from object or V1 -",
	paste0("V", length(covlist)), "\n"))
  }
  if(onames & vnames){
    warning(cat("Variables in calc are mix of names from object and V1 -",
	paste0("V", length(covlist)),
      "\n\tconverting to all V1 -", paste0("V", length(covlist)), "\n"))
    onames <- !onames  #<-- will cause next if statement to give correct names
  }
  if(!onames & vnames) names(covlist) <- covlistVnms


  LHSresult <- eval(deriv(RHS, names(covlist)),
                 covlist)
  gradVec <- matrix(as.vector(attr(LHSresult, "gradient")), ncol=1)
  LHSse <- as.vector(sqrt(crossprod(gradVec, invAI) %*% gradVec))
  LHSresultNm <- if(length(calc) == 3) calc[[2]] else deparse(RHS)

  varcompFunSummary <- data.frame(Estimate = LHSresult, SE = LHSse)
    dimnames(varcompFunSummary) <- list(LHSresultNm, c("Estimate", "Std. Error"))
  class(varcompFunSummary) <- c("deltaSE", class(varcompFunSummary))

 varcompFunSummary
}  #<-- end deltaSE.formula





##############
# List
##############
#' @describeIn deltaSE List method
#' @export
#' @importFrom stats as.formula
deltaSE.list <- function(calc, object, scale = c("theta", "nu")){

  if(!inherits(object, "gremlin")){
    stop(cat("Must supply an object of class", dQuote(gremlin), "\n"))
  }

  cl <- match.call()
  cl_calc <- cl[[match("calc", names(cl))]]
  allFmla <- try(sapply(calc, FUN = inherits, what = "formula"), silent = TRUE)
  if(inherits(allFmla, what = "try-error")){
    warning("calc must be a character expression")
    fmla_lst <- sapply(cl_calc[-1], FUN = function(x){
      as.formula(paste("~", as.expression(x)))})
  } else{
      if(!all(sapply(calc, FUN = inherits, what = "formula"))){
        if(any(sapply(calc, FUN = inherits, what = "formula"))){
          stop("calc must be either all formulas or all expressions")
        }
        fmla_lst <- sapply(calc, FUN = function(x){
          as.formula(paste("~", eval(x)))})
      }

      if(all(sapply(calc, FUN = inherits, what = "formula"))){
        fmla_lst <- calc
      } 
    }  #<-- end if/else error in `allFmla`

    # Last check
    if(any(noForm <- !sapply(fmla_lst, FUN = inherits, what = "formula"))){
      stop(cat("object(s):", which(noForm), "in calc were either not formulas or could not be coerced to formulas\n")) 
    }

  lst_out <- lapply(fmla_lst, FUN = deltaSE.formula, object, scale)

 do.call("rbind", lst_out)
}



