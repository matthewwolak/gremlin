#' Mixed-Effects REML Incorporating Generalized Inverses
#'
#' Fit linear mixed-effects models using restricted (or residual) maximum
#' likelihood (REML) and with generalized inverse matrices to specify covariance
#' structures for random effects. In particular, the package is suited to fit
#' quantitative genetic mixed models, often referred to as 'animal models'
#' (Henderson 1973). Implements the average information algorithm (Johnson &
#' Thompson 1995; Gilmour et al. 1995; Meyer & Smith 1996) as the main tool to
#' maximize the restricted log-likelihood, but with other algorithms available.
#'
#' The average information algorithm combined with sparse matrix techniques can
#' potentially make model fitting very efficient.
#'
#' @aliases gremlin-package
#' @useDynLib gremlin, .registration = TRUE
#' @importFrom methods as is slot
#' @import Matrix
#' @importFrom stats var
#' @references
#'   Henderson, C.R. 1973. Sire evaluation and genetic trends. Journal of Animal
#'   Science 1973:10-41. 
#'   
#'   Johnson, D.L. and R. Thompson. 1995. Restricted maximum likelihood
#'   estimation of variance components for univariate animal models using sparse
#'   matrix techniques and average information. Journal of Dairy Science
#'   78:449-456.
#'
#'   Gilmour, A.R., R. Thompson, and B.R. Cullis. 1995. Average information
#'   REML: An efficient algorithm for variance parameter estimation in linear
#'   mixed models. Biometrics 51:1440-1450.
#'
#'   Meyer, K. and S.P. Smith. 1996. Restricted maximum likelihood estimation for
#'   animal models using derivatives of the likelihood. Genetics, Selection, and
#'   Evolution 28:23-49.
#'
#'   Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding
#'   Values, 2nd ed. CABI Publishing, Cambridge, MA, USA.
#' @examples
#' \dontrun{
#'   # Following the example from Mrode 2005, chapter 11.
#'   library(nadiv)  #<-- to construct inverse of the numerator relatedness matrix
#'   Ainv <- makeAinv(Mrode11[, 1:3])$Ainv
#'
#'   gr11lmm <- gremlin(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	ginverse = list(calf = Ainv),
#'	Gstart = matrix(0.2), Rstart = matrix(0.4),  #<-- specify starting values
#'	maxit = 15,    #<-- maximum iterations
#'      v = 2, vit = 1,  #<-- moderate screen output (`v`) every iteration (`vit`)
#'      algit = "AI")  #<-- only use Average Information algorithm iterations
#'   summary(gr11lmm)
#'
#'   # Compare the model to a Linear Model with no random effects
#'   ## Use `update()` to update the model
#'   gr11lm <- update(gr11lmm, random = ~ 1)  #<-- `~1`=drop all random effects
#'   summary(gr11lm)
#'
#'   # Do analysis of variance between the two models
#'   ## See AIC or evaluate likelihood ratio against a Chi-squared distribution
#'   anova(gr11lm, gr11lmm)
#' }
"_PACKAGE"














 



################################################################################
#' Mixed-effect modeling functions.
#'
#' Fit and setup functions for linear mixed-effect model (Gaussian responses).
#'
#' @aliases gremlin gremlinR gremlinSetup mkModMats update
#' @param formula A \code{formula} for the response variable(s) and fixed effects.
#' @param random A \code{formula} for the random effects.
#' @param rcov A \code{formula} for the residual covariance structure.
#' @param data A \code{data.frame} in which to look for the terms in
#'   \code{formula}, \code{random}, and \code{rcov}.
#' @param subset An expression for the subset of \code{data} to use.
#' @param ginverse A \code{list} of (preferably sparse) inverse matrices that
#'   are proportional to the covariance structure of the random effects.
#'   The name of each element in the list should match a column in \code{data}
#'   that is associated with a random term. All levels of the random term should
#'   appear as \code{rownames} for the matrices.
#' @param Gstart A \code{list} of matrices with starting (co)variance values for
#'   the G-structure or random effects terms.
#' @param Rstart A \code{list} of matrices with starting (co)variance values for
#'   the R-structure or residual terms.
#' @param Bp A prior specification for fixed effects.
#' @param Gcon,Rcon A \code{list} of matrices with constraint codes for the
#'   G-structure/random effects or R-structure/residual effects terms,
#'   respectively. Must be a \code{character} of either \code{"F"} for fixed,
#'   \code{"P"} for positive, or \code{"U"} for unbounded. 
#' @param maxit An \code{integer} specifying the maximum number of likelihood
#'   iterations.
#' @param algit A \code{character} vector of length 1 or more or an expression
#'   to be evaluated that specifies the algorithm to use for proposing
#'   (co)variances in the next likelihood iteration.
#' @param vit An \code{integer} value specifying the verbosity of screen output
#'   on each iteration. A value of zero gives no iteration specific output and
#'   larger values increase the amount of information printed on the screen.
#' @param v An \code{integer} value specifying the verbosity of screen output
#'   regarding the model fitting process. A value of zero gives no details and
#'   larger values increase the amount of information printed on the screen.
#' @param control A \code{list} returned by \code{gremlinControl} containing
#'   specific named values from that function. See \code{\link{gremlinControl}}.
#' @param na.action What to do with NAs.
#' @param offset Should an offset be specified.
#' @param contrasts Specify the type of contrasts for the fixed effects.
#' @param Xsparse Should sparse matrices be used for the fixed effects design
#'   matrix.
#' @param x,object An object of class \code{gremlin}.
#' @param \dots Additional arguments to be passed to control the model fitting.
#'
#' @return A \code{list} containing an object of class \code{grMod} and, if a 
#'   model was fit (\code{gremlin} or \code{gremlinR}) then an object containing
#'   details of the REML iterations (object \code{itMat}). An object of class
#'   \code{grMod} contains:
#'   \describe{
#'     \item{call }{The model \code{call}.}
#'     \item{modMats }{A \code{list} of the model matrices used to construct the
#'       mixed model equations.}
#'       \item{y }{The response vector.}
#'       \item{ny }{The number of responses.}
#'       \item{ncy }{The number of columns of the original response.}
#'       \item{X }{ The fixed effects design matrix.}
#'       \item{nb }{The number of columns in X.}
#'       \item{Zr }{The residual design matrix.}
#'       \item{Zg }{A list of the design matrices for each random term.}
#'       \item{nG }{The number of parameters in the G structure.}
#'       \item{listGeninv }{A list of generalized inverse matrices.}
#'       \item{logDetG }{The log-determinants of the generalized inverse 
#'       matrices - necessary to calculate the log-likelihood.}
#'
#'     \item{rfxIncContrib2loglik }{A \code{numeric} value containing the sum
#'       of the log determinants of the random effects that do not change between
#'       log-likelihood iterations (i.e., the part of the log determinants of 
#'       (co)variance matrices to be estimated that have been factored out).}
#'     \item{ndgeninv }{A \code{logical} indicating which terms in the random
#'       formula have generalized inverses associated with them (non-diagonal 
#'       matrices in the Kronecker product.}
#'     \item{dimsZg, nminffx, rfxlvls, nminfrfx }{\code{Numeric} vectors or scalars
#'       describing the numbers of random effects or some function of random and
#'       fixed effects.}
#'     \item{conv, bounds }{(Co)variance component constraints and boundaries
#'       of the allowable parameter space for each component.}
#'     \item{thetav }{A \code{vector} of the (co)variance parameters to
#'       be estimated by REML withe the attribute \dQuote{skel} giving the
#'       skeleton to recreate a list of \code{matrices} from this vector.}
#'     \item{thetaG, thetaR }{\code{Vectors} indexing the random and residual
#'        (co)variances, respectively, in a list of (co)variance matrices (i.e.,
#'        \code{theta}).}
#'     \item{nu }{A \code{list} of transformed (co)variance matrices
#'       to be fit by REML. If a residual variance has been factored out of the
#'       mixed model equations, \code{nu} contains the \sQuote{lambda}
#'       parameterization with expresses the (co)variance components as ratios
#'       of variance parameters with the residual variance. The \sQuote{nu} scale
#'       (co)variances are the ones actually fit by REML.}
#'     \item{sigma2e }{The estimate of the factored out residual variance from
#'       the mixed model equations (i.e., the \sQuote{lambda} scale)
#'       \eqn{\sigma^{2}_{e}}.}
#'     \item{p }{An \code{integer} for the total number of (co)variances to be
#'       estimated.}
#'     \item{lambda }{A \code{logical} indicating whether the \sQuote{lambda}
#'       scale parameterization has been used.}
#'     \item{uni }{A \code{logical} to indicate if the model is univariate or not.}
#'     \item{W, tWW, RHS, Bpinv }{Sparse matrices of class \code{Matrix} that 
#'       form the mixed model equations and do not change between iterations of
#'       REML. These are the column bound \sQuote{X} and \sQuote{Z} design
#'       matrices for fixed and random effects, the cross-product of \code{W},
#'       the Right-Hand Side of the mixed model equations, and the inverse of
#'       the fixed effect prior matrix (zeroes on the diagonal if no priors have
#'       been specified). Note, these may be \code{NULL} if \code{lambda=FALSE},
#'       because the \code{NULL} objects are not used or do change between
#'       REML iterations.}
#'     \item{sLc }{A \code{Matrix} containing the symbolic Cholesky factorization
#'       of the coefficient matrix of the Mixed Model Equations.}
#'     \item{sln }{A one column \code{matrix} of solutions in the mixed model
#'       equations.}
#'     \item{Cinv_ii }{A one column \code{matrix} of variances for the solutions
#'       to the mixed model equations. These are obtained from the diagonal of
#'       the inverse Coefficient matrix in the mixed model equations. If lambda
#'       is \code{TRUE} then these are on the lambda scale and must be
#'       multiplied by \code{sigma2e} to be converted to the original data scale.}
#'     \item{r }{A one column \code{matrix} of residual deviations, response minus
#'       the values expected based on the solutions, corresponding to the order
#'       in \code{modMats$y}.} 
#'     \item{AI }{A \code{matrix} of values containing the Average Information
#'       matrix, or second partial derivatives of the likelihood with respect to
#'       the transformed (co)variance components (\sQuote{nu}). The inverse of
#'       this matrix gives the sampling (co)variances of these transformed
#'       (co)variance components.}
#'     \item{dLdnu }{A single column \code{matrix} of first derivatives of
#'       the transformed (co)variance parameters (\sQuote{nu}) with respect to
#'       the log-Likelihood.}
#'     \item{maxit }{See the parameter described above.}
#'     \item{algit }{A \code{character} vector of REML algorithms to use in each
#'       iteration.}
#'     \item{vit }{See the parameter described above.}
#'     \item{v }{See the parameter described above.}
#'     \item{cctol }{A \code{numeric} vector of convergence criteria thresholds.
#'       See \code{\link{gremlinControl}} for more details.}
#'     \item{ezero, einf }{\code{numeric} values for the effective numbers to
#'       use as \dQuote{zero} and maximum negative or positive numbers. Values
#'       less than \code{ezero} are treated as zero and fixed to this value.
#'       Values less than \code{-1*einf} or greater than \code{einf} are
#'       restricted to these values. See \code{\link{gremlinControl}} for more
#'       details.}
#'     \item{step }{A \code{numeric} value indicating the step-reduction to use.
#'       See \code{\link{gremlinControl}} for more details.}
#'
#'     \item{itMat }{A \code{matrix} of details about each iteration. Rows
#'       indicate each REML iteration (rownames reflect the REML algorithm used)
#'       and columns contain:
#'       \describe{
#'         \item{nu, theta}{(Co)variance parameters.}
#'         \item{sigma2e }{See \sQuote{sigma2e} described above.}
#'         \item{tyPy, logDetC }{Estimates for two these two components of the
#'           log of the REML likelihoods. These are obtained from Cholesky
#'           factorization of the coefficient matrix of the mixed model equations.}
#'         \item{loglik }{The REML log-likelihood.}
#'         \item{itTime }{Time elapsed for each REML iteration.}
#'       }
#'     }
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#'   grSire <- gremlin(WWG11 ~ sex, random = ~ sire, data = Mrode11)
#'   # Now drop sire random effects and use the `anova` method to compare models
#'   grLM <- update(grSire, random = ~ 1)  #<-- use `~1` to drop all random effects
#`     ## compare models
#'     anova(grSire, grLM)
#'
#'   # Modular functions
#'   ## get model matrices for a mixed model
#'   mM11 <- mkModMats(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#'
#'   ## setup model, but do not evaluate the log-likelihood
#'   grSetup <- gremlinSetup(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#'   ## maximize the restricted maximum likelihood
#'   grOut <- remlIt(grSetup)
#'   summary(grOut)
#'
#' \dontrun{
#'   # Following the example from Mrode 2005, chapter 11.
#'   library(nadiv)  #<-- to construct inverse of the numerator relatedness matrix
#'   Ainv <- makeAinv(Mrode11[, 1:3])$Ainv
#'
#'   gr11lmm <- gremlin(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	ginverse = list(calf = Ainv),
#'	Gstart = matrix(0.2), Rstart = matrix(0.4),  #<-- specify starting values
#'	maxit = 15,    #<-- maximum iterations
#'      v = 2, vit = 1,  #<-- moderate screen output (`v`) every iteration (`vit`)
#'      algit = "AI")  #<-- only use Average Information algorithm iterations
#'   summary(gr11lmm)
#'
#'   # Compare the model to a Linear Model with no random effects
#'   ## Use `update()` to update the model
#'   gr11lm <- update(gr11lmm, random = ~ 1)  #<-- `~1`=drop all random effects
#'   summary(gr11lm)
#'
#'   # Do analysis of variance between the two models
#'   ## See AIC or evaluate likelihood ratio against a Chi-squared distribution
#'   anova(gr11lm, gr11lmm)
#' }
#' @export
gremlin <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL, Bp = NULL,
		Gcon = NULL, Rcon = NULL,
		maxit = 20, algit = NULL,
		vit = 1, v = 1,
		control = gremlinControl(), ...){

  mc <- as.list(match.call())
  mGmc <- as.call(c(quote(gremlinSetup), mc[-1]))
  grMod <- eval(mGmc, parent.frame())
    grMod$call[[1L]] <- quote(gremlin)
    if(v > 2) cat("\tmodel matrices made\n")

  if(v > 2) cat("\tbeginning REML iterations\n")
  grModOut <- remlIt(grMod)

  endTime <- Sys.time()
  if(v > 0) cat("gremlin ended:\t\t", format(endTime, "%H:%M:%S"), "\n")

 return(structure(list(grMod = grModOut$grMod,
		itMat = grModOut$itMat),
	class = c("gremlin"),
	startTime = attr(grMod, "startTime"), endTime = endTime))
}  #<-- end `gremlin()`
################################################################################









#' @rdname gremlin
#' @export
gremlinR <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL, Bp = NULL,
		Gcon = NULL, Rcon = NULL,
		maxit = 20, algit = NULL,
		vit = 1, v = 1,
		control = gremlinControl(), ...){

  mc <- as.list(match.call())
  mGmc <- as.call(c(quote(gremlinSetup), mc[-1]))
  grMod <- eval(mGmc, parent.frame())
    grMod$call[[1L]] <- quote(gremlinR)
  class(grMod) <- c(class(grMod)[1], "gremlinR", class(grMod)[2])
    if(v > 2) cat("\tmodel matrices made\n")

  if(v > 2) cat("\tbeginning REML iterations\n")
  grModOut <- remlIt(grMod)

  endTime <- Sys.time()
  if(v > 0) cat("gremlin ended:\t\t", format(endTime, "%H:%M:%S"), "\n")

 return(structure(list(grMod = grModOut$grMod,
		itMat = grModOut$itMat),
	class = c("gremlinR", "gremlin"),
	startTime = attr(grMod, "startTime"), endTime = endTime))
}  #<-- end `gremlinR()`
################################################################################ 















################################################################################
# UPDATE methods
#################

#' @method getCall gremlin
#' @rdname gremlin
#' @export
getCall.gremlin <- function(x, ...) x$grMod$call
 


###########
# gremlin
###########
#' @method update gremlin
#' @rdname gremlin
#' @export
#' @importFrom stats formula as.formula update.formula
update.gremlin <- function(object, ...){

  call <- getCall(object)
  new_args <- list(...)
  # original `object` class
  oocl <- class(object)

  # First ensure all requested updates to the model are named to prevent chaos
  if(any(names(new_args) == "")){
    stop("All arguments to be updated in the model must be named")
  } 
  # Get names of all defined arguments to `gremlin()` or `gremlinR()`
  defaultCall <- as.list(formals(eval(call[[1L]])))
    # Remove the dots
    dotInd <- match("...", names(defaultCall))
    defaultCall <- defaultCall[-dotInd]   
    grArgNames <- names(defaultCall)
  # Now check all requested updates are named with a valid argument name
  ## first, return index in formal arguments to the function that match
  argNameInd <- pmatch(names(new_args), grArgNames)
  ## make sure all matches worked
  if(any(is.na(argNameInd))){
    noarg <- which(is.na(argNameInd))
    stop(cat("Argument(s)", names(new_args)[noarg], "not valid formal gremlin argument(s)\n"))
  }
  ## replace partial names with full/correct name
  names(new_args) <- grArgNames[argNameInd]
  




  # Go through new_args and determine if a different model is requested
  diffMod <- FALSE
  for(arg in c("formula", "random", "rcov")){
    if((arg %in% names(new_args)) && !is.null(call[[arg]])){
      oldForm <- as.formula(call[[arg]])
      oldFormEnv <- environment(oldForm)
      newForm <- update.formula(oldForm, new_args[[arg]])
      environment(newForm) <- oldFormEnv
  
      if(!identical(oldForm, newForm)){
        diffMod <- TRUE
        if((arg %in% c("random", "rcov")) && (newForm == formula(~1))){
          call[[arg]] <- defaultCall[[arg]]
        } else call[[arg]] <- newForm
      }
    } else{
        if(arg %in% names(new_args) && is.null(call[[arg]])){
          diffMod <- TRUE
          oldForm <- as.formula(defaultCall[[arg]])
          newForm <- update.formula(oldForm, new_args[[arg]])

          if((arg %in% c("random", "rcov")) && (newForm == formula(~1))){
            call[[arg]] <- defaultCall[[arg]]
          } else call[[arg]] <- newForm
        }
      }  #<-- end `else`

  }  #<-- end `for` loop through formulas


  for(arg in c("data", "ginverse", "Bp")){
    if((arg %in% names(new_args)) && !is.null(call[[arg]])){
      if(!identical(eval(call[[arg]]), new_args[[arg]])){
        diffMod <- TRUE
        call[[arg]] <- new_args[[arg]]
      }
    } else{
        if((arg %in% names(new_args)) && is.null(call[[arg]])){
          diffMod <- TRUE
          call[[arg]] <- new_args[[arg]]
        }
      }  #<-- end `else`
  }  #<-- end `for` loop through data/misc objects



  # below don't change the model, just require a re-run with the new values
  if(diffMod){
    if("Gstart" %in% names(new_args)) call[["Gstart"]] <- new_args[["Gstart"]]
    if("Rstart" %in% names(new_args)) call[["Rstart"]] <- new_args[["Rstart"]]
    if("Gcon" %in% names(new_args)) call[["Gcon"]] <- new_args[["Gcon"]]
    if("Rcon" %in% names(new_args)) call[["Rcon"]] <- new_args[["Rcon"]]
  } else{
      # Gstart / Rstart
      ## use output from `object` and replace if G/Rstart provided in `update()`
      theta <- vech2matlist(object$grMod$thetav, attr(object$grMod$thetav, "skel"))
      # Gstart first
      if("Gstart" %in% names(new_args)){
        thetaGtmp <- c(G = sapply(new_args[["Gstart"]], FUN = stTrans))
        if(length(thetaGtmp) != length(object$grMod$thetaG)){
          stop(cat("Gstart of length", length(thetaGtmp),
                   "is not same length as in original model\n"))
        }
        theta[object$grMod$thetaG] <- thetaGtmp  
      } 
      # now Rstart
      if("Rstart" %in% names(new_args)){
        thetaRtmp <- c(R. = stTrans(new_args[["Rstart"]]))
        if(length(thetaRtmp) != length(object$grMod$thetaR)){
          stop(cat("Rstart of length", length(thetaRtmp),
                   "is not same length as in original model\n"))
        }
        theta[object$grMod$thetaR] <- thetaRtmp  
      } 
    
      object$grMod$thetav <- matlist2vech(theta)
      #object$grMod$nu <- theta2nu_trans(theta)
      if(object$grMod$lambda){
        object$grMod$nu <- theta2nu_lambda(theta, object$grMod$thetaG,
                                           object$grMod$thetaR)
      } else object$grMod$nu <- theta

      # Gcon / Rcon
      ## use output from `object` and replace if G/Rcon provided in `update()`
      conv <- object$grMod$conv
      # Gcon first
      if("Gcon" %in% names(new_args)){
        tmpGcon <- conTrans(new_args[["Gcon"]], "P")  #<-- "P"=Rcon placeholder

        if((length(tmpGcon)-1) != length(object$grMod$thetaG)){
          stop(cat("Vectorized Gcon of length", length(tmpGcon)-1,
                   "is not same length as in original model\n"))
        }
        conv[object$grMod$thetaG] <- tmpGcon[object$grMod$thetaG]  
      } 
      # now Rcon
      if("Rcon" %in% names(new_args)){
        tmpRcon <- conTrans("P", new_args[["Rcon"]])  #<-- "P"=Gcon placeholder
        if((length(thetaRtmp)-1) != length(object$grMod$thetaR)){
          stop(cat("Vectorized Rcon of length", length(tmpRcon),
                   "is not same length as in original model\n"))
        }
        conv[object$grMod$thetaR] <- tmpRcon[-1] 
      } 
    
      if(length(object$grMod$conv) != length(conv)){
        stop("Updated constraint vector not same length as 'conv'")
      }
      object$grMod$conv <- conv
      #TODO add levels as add constaints
      # F=Fixed | P=Positive | U=Unbounded | B=Boundary/Bounded
      boundChoices <- matrix(c(NA, NA,                         # Fixed
	             object$grMod$ezero, object$grMod$einf,    # Positive
                     -1*object$grMod$einf, object$grMod$einf), # Unbounded
        ncol = 2, byrow = TRUE)  #<-- TODO add boundaries as add constraints
      object$grMod$bounds[] <- matrix(boundChoices[as.integer(conv), ], ncol = 2)
    
    }  # end if/else diffMod




  # MUST do `maxit` before `algit` to handle `algit` length correctly
  for(arg in c("maxit", "vit", "v")){
    if(arg %in% names(new_args)){
      if(diffMod) call[[arg]] <- new_args[[arg]]
        else object$grMod[[arg]] <- new_args[[arg]]
    }
  }  #<-- end `for` loop through iteration control arguments

  # Handle `algit` separately (must be done after `maxit` and treated special)
  maxitTmp <- c(new_args[["maxit"]], object$grMod$maxit)[which(!sapply(c(new_args[["maxit"]], object$grMod$maxit), FUN = is.null))[1]]
  ## Check validity of `algit`
  ## Fill-in `algit` default if necessary
  if(is.null(new_args[["algit"]])){
    if(is.null(call[["algit"]])) algit <- defaultCall[["algit"]]
      else algit <- call[["algit"]]
    if(is.null(algit)) algit <- c(rep("EM", min(maxitTmp, 2)),
                                  rep("AI", max(0, maxitTmp-2)))
  } else{
      algChoices <- c("EM", "AI", "bobyqa", "NR") #TODO Update if add/subtract any
      algMatch <- pmatch(new_args[["algit"]], algChoices,
        nomatch = 0, duplicates.ok = TRUE)
      if(any(algMatch == 0)){  
        stop(cat("Algorithms:", new_args[["algit"]][which(algMatch == 0)],
        "not valid. Please check values given to the `algit` argument\n"))
      }
      algit <- algChoices[algMatch]
    }  #<-- end if/else new_args null for algit
  if(length(algit) == 1) algit <- rep(algit, maxitTmp)
  if(length(algit) > maxitTmp) algit <- algit[1:maxitTmp]
  if(length(algit) < maxitTmp) algit <- rep(tail(algit, 1), maxitTmp)
  if(diffMod) call[["algit"]] <- algit
    else object$grMod[["algit"]] <- algit


  # gremlinControl() changes
  ##TODO add `algorithm` and `algArgs`?
  if("control" %in% names(new_args)){
    if(diffMod){
      call[["control"]] <- new_args[["control"]]
    } else{
        for(arg in c("cctol", "ezero", "einf", "step", "lambda")){
          object$grMod[arg] <- new_args$control[arg]
        }  #<-- end for arg
      }  #<-- end if/else diffMod
  }  #<-- end if


  ################################################################


  if(diffMod){
    mGmc <- as.call(c(quote(gremlinSetup), as.list(call)[-1]))
    grMod <- eval(mGmc, parent.frame())
      if(is(oocl, "gremlinR")){
        grMod$call[[1L]] <- quote(gremlinR)
        class(grMod) <- c(class(grMod)[1], "gremlinR", class(grMod)[2])
      } else grMod$call[[1L]] <- quote(gremlin)

    grModOut <- remlIt(grMod)

  } else{
      # need to re-set `sln` and `Cinv_ii` to zeroes so c++ calculations work
      object$grMod$sln[] <- object$grMod$Cinv_ii[] <- rep(0, length(object$grMod$sln))
      grModOut <- remlIt(object$grMod)
        # If model has not changed PREpend previous itMat to latest
        grModOut$itMat <- rbind(object$itMat, grModOut$itMat)
        #TODO not sure what to do about attr(*, "startTime"): 
        ##Which to use? Or include both somehow? Now defaults to start of update

    }  #<-- end if/else diffMod



  endTime <- Sys.time()
  if(grModOut$grMod$v > 0) cat("gremlin ended:\t\t", format(endTime, "%H:%M:%S"), "\n")

 return(structure(list(grMod = grModOut$grMod,
		itMat = grModOut$itMat),
	class = oocl,
	startTime = attr(grModOut$grMod, "startTime"), endTime = endTime))
}  #<-- end `update.gremlin`











################################################################################
#' @rdname gremlin
#' @export
#' @import Matrix
gremlinSetup <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL, Bp = NULL,
		Gcon = NULL, Rcon = NULL,
		maxit = 20, algit = NULL,
		vit = 1, v = 1,
		control = gremlinControl(), ...){

  stopifnot({
    inherits(formula, "formula")
    length(formula) == 3L
    inherits(control, "gremlinControl")
  })
  mc <- as.list(match.call())
  startTime <- Sys.time()
  if(v > 0) cat("gremlin started:\t\t", format(startTime, "%H:%M:%S"), "\n")
  m <- match(c("formula", "random", "rcov", "data", "subset", "ginverse", "na.action", "offset", "contrasts", "Xsparse"), names(mc), 0)
  mMmc <- as.call(c(quote(mkModMats), mc[m]))
  modMats <- eval(mMmc, parent.frame())

  if(missing(rcov)) mc$rcov <- as.list(formals(eval(mc[[1L]])))[["rcov"]]
  #algChoices <- c("EM", "AI", "bobyqa", "NR", control$algorithm)
  algChoices <- c("EM", "AI", "bobyqa", "NR")  #<-- ignore control$algorithm
    if(!is.null(control$algorithm)){
      #TODO check validity of `control$algorithm` and `control$algArgs`
      ## need to pass algorithm to `gremlinR` or switch to it if `gremlin` called
      ## temporarily IGNORE with warning
      warning(cat("Ignored algorithm(s) supplied in control",
        dQuote(control$algorithm),
        ". gremlin is not old enough for user-specified algorithms\n"),
          immediate. = TRUE)
    }
  algMatch <- pmatch(algit, algChoices, nomatch = 0, duplicates.ok = TRUE)
  if(all(algMatch == 0) & !all(algit %in% control$algorithm)){ 
      stop(cat("Algorithms:", dQuote(algit[which(algMatch == 0)]),
        "not valid. Please check values given to the `algit` argument\n"))
  }
  if(any(algMatch == 0)){  
    warning(cat("Algorithms:", dQuote(algit[which(algMatch == 0)]),
      "not valid - dropped from the list\n"))
    algit <- algit[-which(algMatch == 0)]
  }
  if(is.null(mc$algit)){
    algit <- c(rep("EM", min(maxit, 2)), rep("AI", max(0, maxit-2)))
  } else algit <- algChoices[algMatch]
  if(length(algit) == 0) algit <- c(rep("EM", min(maxit, 2)),
                                    rep("AI", max(0, maxit-2)))
  if(length(algit) == 1) algit <- rep(algit, maxit)



  #TODO check dimensions G/Rstart
#FIXME assumes univariate
  if(modMats$nG > 8) fra <- 0.9 / modMats$nG else fra <- 0.1
  if(is.null(mc$Gstart)) Gstart <- as.list(rep(fra*var(modMats$y), modMats$nG))
    else Gstart <- eval(mc$Gstart)
#FIXME assumes univariate
  if(is.null(mc$Rstart)) Rstart <- matrix(0.1*var(modMats$y))
    else Rstart <- eval(mc$Rstart)
#FIXME assumes univariate
  if(is.null(mc$Gcon)) Gcon <- as.list(rep("P", modMats$nG))
    else Gcon <- eval(mc$Gcon)
#FIXME assumes univariate
  if(is.null(mc$Rcon)) Rcon <- matrix("P")
    else Rcon <- eval(mc$Rcon)

  thetaSt <- start2theta(Gstart, Rstart, name = names(modMats$Zg))
  thetav <- matlist2vech(thetaSt$theta)
  p <- length(thetav)
  conv <- conTrans(Gcon, Rcon)
    if(length(conv) != p){
      if(length(Gcon) != length(Gstart)){
        stop(cat("Gcon (length=", length(Gcon),
          ") must be same length as Gstart (", length(Gstart), ")\n"))
      }
      if(length(Rcon) != length(Rstart)){
        stop(cat("Rcon (length=", length(Rcon),
          ") must be same length as Rstart (", length(Rstart), ")\n"))
      }

      GconLengths <- sapply(Gcon, FUN = length)
      GstartLengths <- sapply(Gstart, FUN = length)
      RconLengths <- sapply(Rcon, FUN = length)
      RstartLengths <- sapply(Rstart, FUN = length)
      if(any((GconLengths - GstartLengths) != 0)){
        stop(cat("Gcon element(s)", which((GconLengths - GstartLengths) != 0),
          "are not of the same length as their corresponding Gstart elements\n"))
      }
      if(any((RconLengths - RstartLengths) != 0)){
        stop(cat("Rcon element(s)", which((RconLengths - RstartLengths) != 0),
          "are not of the same length as their corresponding Rstart elements\n"))
      }
    }
    names(conv) <- names(thetav)
    #TODO add levels as add constaints
    # F=Fixed | P=Positive | U=Unbounded | B=Boundary/Bounded
    conv <- factor(conv, levels = c("F", "P", "U", "B"))
  boundChoices <- matrix(c(NA, NA,                         # Fixed
			   control$ezero, control$einf,    # Positive
                           -1*control$einf, control$einf), # Unbounded
      ncol = 2, byrow = TRUE)  #<-- TODO add boundaries as add constraints
  bounds <- matrix(boundChoices[as.integer(conv), ], ncol = 2,
    dimnames = list(names(conv), c("LB", "UB")))









#TODO TODO TODO TODO TODO TODO TODO
#XXX Determine from model how to make `lambda==FALSE` for models where appropriate
  lambda <- control$lambda







#TODO put `uni` in `mkModMats()`
  if(modMats$ncy == 1) uni <- TRUE
    else stop("gremlin isn't old enough to play with multivariate models")




#FIXME: change G to cholesky of G with log(diagonals)
## e.g., parameterisation to ensure postive-definiteness
### when to do cholesky factorization of G?
#### is it easier to do direct product of cholesky's than direct product then cholesky?
#### If the latter then save the symbolic factorization and just do updates?
#TODO: when change to cholesky with log(diags) - change also in `update.gremlin()`
  #nu <- theta2nu_trans(thetaSt$theta)
  if(lambda) nu <- theta2nu_lambda(thetaSt$theta, thetaSt$thetaG, thetaSt$thetaR)
    else nu <- thetaSt$theta


    # 1 Create coefficient matrix of MME (C)
    ##1a form W by `cbind()` X and each Z_i
    W <- modMats$X
    for(g in seq_len(modMats$nG)) W <- cbind(W, modMats$Zg[[g]])

    # Rand Fx incidence matrix part of 'log(|G|)'
    #FIXME: Only works for independent random effects right now!
    rfxIncContrib2loglik <- sum(unlist(modMats$logDetG))
    if(is.null(Bp)){
      Bp <- as(diag(x = 0, nrow = modMats$nb, ncol = modMats$nb), "dgCMatrix")
    } else{
       #TODO check and maybe create prior from Bp specified in call
       stop("Currently can't take fixed effect prior")
      } 
    if(all(Bp@x == 0)){
      Bpinv <- Bp + diag(0, nrow(Bp))
        if(length(Bpinv@x) == 0){
          # need to put explicit 0s on diagonal
          Bpinv@x <- as.double(rep(0, nrow(Bp)))
          Bpinv@i <- as.integer(seq(nrow(Bp))-1)
          Bpinv@p <- as.integer(c(seq(nrow(Bp))-1, nrow(Bp)))
        }
    } else Bpinv <- solve(Bp)
      # Bpinv <-- used every iteration
      ## `Bpinv` replaces `zero` in earlier version of `gremlinR()`
      ## Can allow for prior on fixed effects
      ## (see Schaeffer 1991 summary of Henderson's result)
    # if no random effects in model:
    if(modMats$nG < 1){
      ndgeninv <- c(0)
      dimsZg <- matrix(0, nrow = 2, ncol = 1)
    } else{
        # Find Non-diagonal ginverses
        ## FALSE/0=I; TRUE/1=geninv 
        ndgeninv <- sapply(seq_len(modMats$nG),
	  FUN = function(g){class(modMats$listGeninv[[g]]) != "ddiMatrix"}) 
        dimsZg <- sapply(seq_len(modMats$nG),
	  FUN = function(g){slot(modMats$Zg[[g]], "Dim")})
      }
    sln <- Cinv_ii <- matrix(0, nrow = modMats$nb + sum(dimsZg[2, ]), ncol = 1)
    r <- matrix(0, nrow = modMats$ny, ncol = 1)
    if(lambda){
      tWW <- crossprod(W)
      RHS <- Matrix(crossprod(W, modMats$y), sparse = TRUE)  # <-- Same every iteration
    } else tWW <- RHS <- NULL
    sLc <- NULL  #<-- initialize NULL and will generate vs. update if is.null()

#TODO put these with `mkModMats()` - need to figure out multivariate version/format
    # 5b log(|R|) and log(|G|) <-- Meyer 1989 (uni) & 1991 (multivar)
    # Only have to do these once per model
#FIXME make sure `nminffx` == `ncol(X)` even when reduced rank
    nminffx <- modMats$ny - modMats$nb
    rfxlvls <- sapply(modMats$Zg, FUN = ncol)  #<-- =`dimsZg[2,]` but this gives names
    nr <- if(length(rfxlvls) == 0) 0 else sum(rfxlvls)
    nminfrfx <- nminffx - nr

    AI <- matrix(NA, nrow = p, ncol = p)
    dLdnu <- matrix(NA, nrow = nrow(AI), ncol = 1,
      dimnames = list(names(thetav), NULL))
    sigma2e <- if(lambda) numeric(1) else NA #<-- only needed for `lambda` model

 return(structure(list(call = as.call(mc),
		modMats = modMats,
		rfxIncContrib2loglik = rfxIncContrib2loglik,
		ndgeninv = ndgeninv, dimsZg = dimsZg, nminffx = nminffx,
		rfxlvls = rfxlvls, nminfrfx = nminfrfx,
		conv = conv, bounds = bounds,
		thetav = thetav,
		thetaG = thetaSt$thetaG, thetaR = thetaSt$thetaR,
		nu = nu,
		sigma2e = sigma2e,
		p = p, lambda = lambda, uni = uni,
		W = W, tWW = tWW, RHS = RHS, Bpinv = Bpinv,
		sLc = sLc,
		sln = sln, Cinv_ii = Cinv_ii, r = r,
		AI = AI, dLdnu = dLdnu,
		maxit = maxit, algit = algit, vit = vit, v = v,
		cctol = control$cctol,
		ezero = control$ezero, einf = control$einf, step = control$step),
	class = c("grMod", "gremlin"),
	startTime = startTime))
}  #<-- end `gremlinSetup()`
################################################################################
























################################################################################
#' @method is gremlin
#' @rdname gremlin
#' @export
is.gremlin <- function(x) inherits(x, "gremlin")

#' @method is grMod
#' @rdname gremlin
#' @export
is.grMod <- function(x) inherits(x, "grMod")





 








