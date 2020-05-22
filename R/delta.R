################################################################################
# Delta Method methods
#######################

##############
# Generic
##############
#' Standard Error Calculation for Functions of (co)variance parameters.
#'
#' Calculates the standard error for results of simple mathematical functions of
#'   (co)variance parameters using the delta method (Lynch and Walsh 1998,
#'   Appendix 1).
#'
#' 
#'
#' @aliases deltaSE deltaSE.default deltaSE.formula deltaSE.list
#' @param fmla,expr,lst An \code{expression}, \code{formula}, or list (of
#'   \code{formula} or \code{expression}) expressing the mathematical function
#'   of (co)variance component for which to calculate standard errors.
#' @param object A fitted model object of \code{class} \sQuote{gremlin}.
#' @param scale A \code{character} indicating whether to calculate the function
#'   and standard error on the original data scale (\code{\dQuote{theta}}) or
#'   on the underlying scale to which (co)variance components are transformed
#'   to which the model fitting calculations occur (\code{dQuote{nu}}). Defaults
#'   to \code{\dQuote{theta}} if not specified.
#' @return A \code{data.frame} containing the \dQuote{Estimate} and
#'   \dQuote{Std. Error} for the mathematical function of (co)variance components.
#' @references 
#'   TODO (Lynch and Walsh 1998)
#' @author \email{matthewwolak@@gmail.com}
#' @seealso First derivatives are calculated by symbolic differentiation with
#'   \code{\link{deriv}}
#' @import stats:::deriv
#' @examples
#'   # Calculate the sum of the variance components 
#'     grS <- gremlin(WWG11 ~ sex, random = ~ sire, data = Mrode11)
#'     deltaSE(Vsum ~ V1 + V2, grS)
#'
#'   # Calculate standard deviations (with standard errors) from variances
#'     ## Uses a list argument to `fmla`
#'     deltaSE(c(SD1 ~ sqrt(V1), SDresid ~ sqrt(V2)), grS)
#'     deltaSE(list("sqrt(V1)", "sqrt(V2)"), grS)
#'     deltaSE(list(sqrt(V1), sqrt(V2)), grS)
#'
#'   # Additive Genetic Variance calculated from observed Sire Variance
#'     ## First simulate Full-sib data
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
#'     dOut <- deltaSE.formula(thetaV1 ~ V1*V2, grSl, "nu")  #<-- V2 is sigma2e
#'     aiFnOut <- nuVar2thetaVar_lambda(grSl)[1]  #<-- variance (do sqrt below)
#'     stopifnot(abs(deltaOut[, "Std. Error"] - sqrt(covFunOut)) < 1e-10)
#'
#' @export
deltaSE <- function(expr, object, scale = c("theta", "nu")){
  UseMethod("deltaSE", expr)
}

##############
# Default
##############
#' @describeIn deltaSE Default method
deltaSE.default <- function(expr, object, scale = c("theta", "nu")){
  if(!inherits(object, "gremlin")){
    stop(cat("Must supply an object of class", dQuote(gremlin), "\n"))
  }

  fmla <- as.formula(paste("~", as.expression(eval(expr, parent.frame()))))

 deltaSE.formula(fmla, object, scale)
}  #<-- end deltaSE.default





deltaSE.list(list(SD1 ~ sqrt(V1), SD2 ~ sqrt(V2)), grS)
deltaSE.list(list("sqrt(V1)", "sqrt(V2)"), grS)
deltaSE.list(list(sqrt(V1), sqrt(V2)), grS)




##############
# Formula
##############
#' @describeIn deltaSE Formula method
deltaSE.formula <- function(fmla, object, scale = c("theta", "nu")){
  if(!inherits(object, "gremlin")){
    stop(cat("Must supply an object of class", dQuote(gremlin), "\n"))
  }

  RHS <- fmla[[length(fmla)]]

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
  onames <- any(names(covlist) %in% as.list(RHS))
  vnames <- any(paste0("V", seq(length(covlist))) %in% as.list(RHS))
  if(!onames && !vnames){
    stop("Variables in expr/fmla not names from object or V1...Vn")
  }
  if(onames & vnames){
    warning(cat("Variables in expr/fmla are mix of names from object and V1...Vn",
      "\n\tconverting to all V1...Vn"))
    onames <- !onames  #<-- will cause next if statement to give correct names
  }
  if(!onames & vnames) names(covlist) <- paste0("V", seq(length(covlist)))

  LHSresult <- eval(stats:::deriv(RHS, names(covlist)),
                 covlist)
  gradVec <- matrix(as.vector(attr(LHSresult, "gradient")), ncol=1)
  LHSse <- as.vector(sqrt(crossprod(gradVec, invAI) %*% gradVec))
  LHSresultNm <- if(length(fmla) == 3) fmla[[2]] else deparse(RHS)

  varcompFunSummary <- data.frame(Estimate = LHSresult, SE = LHSse)
    dimnames(varcompFunSummary) <- list(LHSresultNm, c("Estimate", "Std. Error"))
  class(varcompFunSummary) <- c("deltaSE", class(varcompFunSummary))

 varcompFunSummary
}  #<-- end deltaSE.formula


##############
# List
##############
#' @describeIn deltaSE List method
deltaSE.list <- function(lst, object, scale = c("theta", "nu")){
  if(!inherits(object, "gremlin")){
    stop(cat("Must supply an object of class", dQuote(gremlin), "\n"))
  }

  cl <- match.call()
  cl_lst <- cl[[match("lst", names(cl))]]
  allFmla <- try(sapply(lst, FUN = inherits, what = "formula"), silent = TRUE)
  if(inherits(allFmla, what = "try-error")){
    fmla_lst <- sapply(cl_lst[-1], FUN = function(x){
      as.formula(paste("~", as.expression(x)))})
  } else{
      if(!all(sapply(lst, FUN = inherits, what = "formula"))){
        if(any(sapply(lst, FUN = inherits, what = "formula"))){
          stop("lst must be either all formulas or all expressions")
        }
        fmla_lst <- sapply(lst, FUN = function(x){
          as.formula(paste("~", as.expression(eval(x, parent.frame()))))})
      }

      if(all(sapply(lst, FUN = inherits, what = "formula"))){
        fmla_lst <- lst
      } 
    }  #<-- end if/else error in `allFmla`

    # Last check
    if(any(noForm <- !sapply(fmla_lst, FUN = inherits, what = "formula"))){
      stop(cat("object(s):", which(noForm), "in lst were either not formulas or could not be coerced to formulas\n")) 
    }

  lst_out <- lapply(fmla_lst, FUN = deltaSE.formula, object, scale)

 do.call("rbind", lst_out)
}


deltaSE.list(list(SD1 ~ sqrt(V1), SD2 ~ sqrt(V2)), grS)
deltaSE.list(list("sqrt(V1)", "sqrt(V2)"), grS)
deltaSE.list(list(sqrt(V1), sqrt(V2)), grS)









