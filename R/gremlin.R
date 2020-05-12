#' Mixed-Effects REML Incorporating Generalized Inverses
#'
#' Fit linear mixed-effects models using restricted (or residual) maximum
#' likelihood (REML) and with generalized inverse matrices to specify covariance
#' structures for random effects. In particular, the package is suited to fit
#' quantitative genetic mixed models, often referred to as 'animal models'.
#' Implements the average information algorithm as the main tool to maximize the
#' restricted likelihood, but with other algorithms available.
#'
#' The package also implements the average information algorithm to efficiently
#' maximize the log-likelihood (Thompson & Johnson 1995; Gilmour et al. 1995;
#' Meyer & Smith 1996). The average information algorithm combined with sparse
#' matrix techniques can potentially make model fitting very efficient.
#'
#' @aliases gremlin-package
#' @useDynLib gremlin, .registration = TRUE
#' @importFrom methods as is slot
#' @import Matrix
#' @importFrom stats var
#' @references
#'   Mrode. 2005.
#'   Meyer & Smith. 1996.
#'   Gilmour et al. 1995.
#'   Thompson & Johnson. 1995.
#' @seealso \code{\link[MCMCglmm:MCMCglmm-package]{MCMCglmm}}
#' @examples
#' \dontrun{
#'   # Following the example from Mrode 2005, chapter 11.
#'   library(nadiv)  #<-- to construct inverse of the numerator relatedness matrix
#'   Ainv <- makeAinv(Mrode11[, 1:3])$Ainv
#'
#'   grOut11 <- gremlin(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	ginverse = list(calf = Ainv),
#'	Gstart = matrix(0.2), Rstart = matrix(0.4),  #<-- specify starting values
#'	maxit = 10,    #<-- maximum iterations
#'      v = 2, vit = 1,  #<-- moderate screen output (`v`) every iteration (`vit`)
#'      algit = "AI")  #<-- only use Average Information algorithm iterations
#'   summary(grOut11)
#'
#'   # Compare the model to a Linear Model with no random effects
#'   ## Use `update()` to update the model
#'   grOutLM <- update(grOut11, random = ~ 1)  #<-- `~1`=drop all random effects
#'   summary(grOutLM)
#'
#'   # Do analysis of variance between the two models
#'   ## See AIC or evaluate likelihood ratio against a Chi-squared distribution
#'   anova(grOutLM, grOut11)
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
#' @param Gstart A \code{list} of starting (co)variance values for the
#'   G-structure or random effects terms.
#' @param Rstart A \code{list} of starting (co)variance values for the
#'   R-structure or residual terms.
#' @param Bp A prior specification for fixed effects.
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
#' @param na.action What to do with NAs.
#' @param offset Should an offset be specified.
#' @param contrasts Specify the type of contrasts for the fixed effects.
#' @param Xsparse Should sparse matrices be used for the fixed effects design
#'   matrix.
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
#'     \item{thetav }{A \code{vector} of the (co)variance parameters to
#'       be estimated by REML withe the attribute \dQuote{skel} giving the
#'       skeleton for recreating a list of \code{matrices} from this vector.}
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
#'       REML. These are the column binded \sQuote{X} and \sQuote{Z} design
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
#'       to the mixed model equations. These are obtained from the diagonal of #'       the inverse Coefficient matrix in the mixed model equations. If lambda
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
#'     \item{cctol }{A \code{numeric} vector of convergence criteria thresholds.}
#'     \item{ezero }{A \code{numeric} value for the effective number to use as
#'       \dQuote{zero}. VAlues less than this number are treated as zero and
#'       fixed to this value.}
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
#' @references
#' Meyer.
#' Henderson.
#' Mrode. 2005.
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
#' @export
gremlin <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL, Bp = NULL,
		maxit = 20, algit = NULL,
		vit = 10, v = 1,
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
		maxit = 20, algit = NULL,
		vit = 10, v = 1,
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
#' @rdname gremlin
#' @export
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
  for(arg in c("Gstart", "Rstart")){
    if(arg %in% names(new_args)){
      if(diffMod) call[[arg]] <- new_args[[arg]]
        else{
          ## fill in G/Rstart if arg in new_args or if not in new_args then fill object starts/thetav/nu with last parameter values to start from there? (Find out what remlIt needs)
          #TODO G/Rstarts need to be transformed to replace thetav/nu, maybe sigma2e, etc. in `object$grMod`

        }  #<-- end `else`
    }  #<-- end if arg %in% names
  }  #<-- end `for` loop through G/Rstart


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
		maxit = 20, algit = NULL,
		vit = 10, v = 1,
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
  if(is.null(mc$Gstart)) Gstart <- as.list(rep(0.1*var(modMats$y), modMats$nG))
    else Gstart <- eval(mc$Gstart)
#FIXME assumes univariate
  if(is.null(mc$Rstart)) Rstart <- matrix(0.5*var(modMats$y))
    else Rstart <- eval(mc$Rstart)

  thetaSt <- start2theta(Gstart, Rstart, name = names(modMats$Zg))
  thetav <- matlist2vech(thetaSt$theta)
  p <- length(thetav)



#TODO TODO TODO TODO TODO TODO TODO
#XXX Determine from model how to make `lambda==FALSE` for models where appropriate
#TODO make `lambda == FALSE` if EM is part of algit (not sure if EM works for lambda or not) #<-- FIXME figure out how to do EM on lambda
  lambda <- control$lambda


#TODO put `uni` in `mkModMats()`
  if(modMats$ncy == 1) uni <- TRUE
    else stop("gremlin isn't old enough to play with multivariate models")




#FIXME: change G to cholesky of G with log(diagonals)
## e.g., parameterisation to ensure postive-definiteness
### when to do cholesky factorization of G?
#### is it easier to do direct product of cholesky's than direct product then cholesky?
#### If the latter then save the symbolic factorization and just do updates?
#XXX Don't need for EM algorithm
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
		cctol = control$cctol, ezero = control$ezero,
		step = control$step),
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





 











################################################################################
################################################################################
###########
# Generic
###########
################################################################################
#' Mixed-effect model Restricted Maximum Likelihood (REML) iterations.
#'
#' Conduct REML iterations to estimate (co)variance parameters of a linear
#'   mixed-effect model (Gaussian responses).
#'
#' @aliases remlIt remlIt.default remlIt.gremlinR
#' @param grMod A gremlin model of class \code{grMod}. See \code{\link{gremlin}}
#'   or \code{\link{gremlinSetup}} for the functions constructing an object
#'   of class \code{grMod}.
#' @param \dots Additional arguments to be passed to control the model fitting.
#'
#' @return A \code{list} containing an object of class \code{grMod} and
#'   \code{matrix} containing details of the REML iterations (object
#'   \code{itMat}). See \code{\link{gremlin}} for descriptions of \code{grMod}
#'   and \code{itMat} objects.
#'
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#'  TODO XXX XXX Make example of each modular component working with next
#'   mod11 <- gremlinSetup(WWG11 ~ sex - 1,
#'   	random = ~ calf,
#'   	data = Mrode11,
#'   	Gstart = matrix(0.2), Rstart = matrix(0.4),
#'   	maxit = 10, v = 2)
#'   REMLout <- remlIt(mod11)
#'
#' @export
remlIt <- function(grMod, ...){
  UseMethod("remlIt", grMod)
}




###########
# Default
###########
#' @describeIn remlIt Default method
#' @export
#' @import Matrix
remlIt.default <- function(grMod, ...){

  gnu <- lapply(grMod$nu, FUN = as, "dgCMatrix") #FIXME do this directly to begin with or just use dense matrices (class="matrix")
  nuv <- matlist2vech(grMod$nu)
  # convert algorithms for each iteration into integers
  intfacalgit <- as.integer(factor(grMod$algit[1:grMod$maxit],
    levels = c("EM", "AI"), ordered = TRUE))
    ## Check that all are implemented in cpp (currently only AI or EM)
    if(any(is.na(intfacalgit))){
      stop(cat("Algorithm", grMod$algit[which(is.na(intfacalgit))],
	"not implemented in c++, try", sQuote(`gremlinR()`), "\n"))
    }

  nnzWG <- with(grMod, c(length(W@x),		# No. nonzero W
    sapply(seq_len(length(thetaG)),
	    FUN = function(g){length(modMats$listGeninv[[g]]@x)}))) # No. nz geninvs

  dimZWG <- with(grMod, c(rowSums(dimsZg),		# Z Dims
    W@Dim))						# W Dims
  if(!is.null(grMod$modMats$listGeninv)){  		
    dimZWG <- c(dimZWG, unlist(lapply(grMod$modMats$listGeninv,  # geninv Dims
      FUN = function(g){g@Dim[[1L]]})))
    geninv_i <- unlist(lapply(grMod$modMats$listGeninv[grMod$ndgeninv],
	FUN = function(g){g@i}))
    geninv_p <- unlist(lapply(grMod$modMats$listGeninv[grMod$ndgeninv],
	FUN = function(g){g@p}))
    geninv_x <- unlist(lapply(grMod$modMats$listGeninv[grMod$ndgeninv],
	FUN = function(g){g@x}))
  } else{
      dimZWG <- c(dimZWG, 0)
      geninv_i <- geninv_p <- geninv_x <- 0
    warning(cat("\n\ngremlin's c++ code is far too sophisticated for such a simple model.\n",
       "\tRe-fitting with the", sQuote('gremlinR()'), "function instead\n\n"))
      return(remlIt.gremlinR(grMod))
    }

  Cout <- .C("ugremlin", PACKAGE = "gremlin",
	as.double(grMod$modMats$y),
	as.integer(grMod$modMats$ny),
	as.integer(grMod$nminffx),		# No. observations - No. Fxd Fx
	as.integer(grMod$ndgeninv),		# non-diagonal ginverses
	as.integer(c(grMod$dimsZg)),
	as.integer(dimZWG), 			# Z, W, and geninv Dims
	as.integer(nnzWG), 			# No. nnz in W and geninvs
	as.integer(grMod$W@i), 					     #W
	as.integer(grMod$W@p),
	as.double(grMod$W@x),
	as.integer(geninv_i), #geninv (generalized inverses)

	as.integer(geninv_p),
	as.double(geninv_x),
	as.double(grMod$rfxIncContrib2loglik),		# Random Fx contribution to log-Likelihood
        as.integer(grMod$lambda),		# TRUE/FALSE lambda (variance ratio?)
	as.integer(grMod$p),				#p=No. nu params
	as.integer(c(length(grMod$thetaG), length(grMod$thetaR))), #No. G and R nus
	as.integer(unlist(lapply(gnu, FUN = function(g) g@Dim[[1L]]))),#dim GRs
	as.integer(unlist(lapply(gnu, FUN = function(g) g@i))),	      #i GRs
	as.integer(unlist(lapply(gnu, FUN = function(g) g@p))),	      #p GRs
	as.integer(unlist(lapply(gnu, FUN = function(g) length(g@x)))), #no. non-zero GRs
	as.double(unlist(lapply(gnu, FUN = function(g) g@x))),	#nu vector
	as.integer(length(grMod$Bpinv@x)),		#Bpinv (fixed fx prior inverse)
	as.integer(grMod$Bpinv@i),
	as.integer(grMod$Bpinv@p),
	as.double(grMod$Bpinv@x),
	as.double(rep(0, grMod$p)),			#empty dLdnu
	as.double(rep(0, grMod$p^2)),	#empty column-wise vector of AI matrix
	as.double(c(grMod$sln)), 			#empty sln
        as.double(c(grMod$Cinv_ii)),			#empty diag(Cinv)
	as.double(c(grMod$r)),				#empty resdiuals
	as.double(rep(0, grMod$maxit*(grMod$p+5))),	#itMat
	as.integer(intfacalgit -1), 			#algorithm for each iteration
	as.integer(grMod$maxit),			#max it./n algit
	as.double(grMod$step),			#init./default step-halving value
	as.double(grMod$cctol),				#convergence tol.
	as.double(grMod$ezero),				#effective 0
#uni?
	as.integer(grMod$v),				#verbosity
	as.integer(grMod$vit),				#when to output status
	as.integer(rep(0, length(grMod$sln))))		#empty sLc->pinv

  i <- Cout[[34]]  #<-- index from c++ always increments +1 at end of for `i`

  grMod$nu[] <- vech2matlist(Cout[[22]], attr(grMod$thetav, "skel"))
  grMod$dLdnu[] <- Cout[[27]]
  if(all(Cout[[28]] == 0)) grMod$AI <- NULL else{
    grMod$AI <- matrix(Cout[[28]], nrow = grMod$p, ncol = grMod$p, byrow = FALSE)
      dimnames(grMod$AI) <- list(rownames(grMod$dLdnu), rownames(grMod$dLdnu))
  }
  grMod$sln[] <- Cout[[29]]
  grMod$Cinv_ii <- Cout[[30]] 
  grMod$r[] <- Cout[[31]]
  #TODO Will definitely need R vs. c++ methods for `update.gremlin()`
  #### can directly use R's `grMod$sLc`, but will need to figure out how to give c++'s `cs_schol()` a pinv (need to reconstruct `sLc` in c++ around pinv (see old code on how I may have done this when I made sLc from sLm)
  grMod$sLcPinv <- Cout[[39]]

  itMat <- matrix(Cout[[32]][1:(i*(grMod$p+5))], nrow = i, ncol = grMod$p+5,
           byrow = TRUE)
    dimnames(itMat) <- list(paste(seq(i), c("EM", "AI")[Cout[[33]][1:i] + 1],
                sep = "-"),
	    c(paste0(names(nuv), "_nu"), "sigma2e",
               "tyPy", "logDetC", "loglik", "itTime"))

  if(grMod$lambda){
    itMat <- cbind(itMat,
      t(apply(itMat[, c(paste0(names(nuv), "_nu"), "sigma2e"), drop = FALSE],
	MARGIN = 1,
	FUN = function(itvec){ matlist2vech(nu2theta_lambda(itvec[1:grMod$p],
          sigma2e = itvec[grMod$p+1], grMod$thetaG, grMod$thetaR))})))[, c(seq(grMod$p),
      seq(grMod$p+6, 2*grMod$p+5), seq(grMod$p+1, grMod$p+5)), drop = FALSE] #<-- thetas named 'nu' in colnames for now, use numeric indices to rearrange

  } else{
      itMat <- cbind(itMat,
        t(apply(itMat[, c(paste0(names(nuv), "_nu"), "sigma2e"), drop = FALSE],
	  MARGIN = 1,
	  FUN = function(itvec){ matlist2vech(nu2theta_noTrans(itvec[1:grMod$p],
	    grMod$thetaG, grMod$thetaR))})))[, c(seq(grMod$p),
        seq(grMod$p+6, 2*grMod$p+5), seq(grMod$p+1, grMod$p+5)), drop = FALSE] #<-- thetas named 'nu' in colnames for now, use numeric indices to rearrange

    }  #<-- end if/else lambda
  # Now sort out the column names
  colnames(itMat) <- c(paste0(names(nuv), "_nu"),
			paste0(names(nuv), "_theta"),
			"sigma2e", "tyPy", "logDetC", "loglik", "itTime")

  grMod$sigma2e[] <- itMat[i, "sigma2e"]
  grMod$thetav[] <- itMat[i, paste0(names(nuv), "_theta")]

 return(structure(list(grMod = grMod,
		itMat = itMat),
	class = "gremlin"))
}  #<-- end `remlIt.default()`
















#################
# R-based method
#################
#' @describeIn remlIt gremlinR method
#' @export
#' @import Matrix
remlIt.gremlinR <- function(grMod, ...){
  # pull a few objects out that will be used repeatedly
  ## favor "small" objects. keep large objects in grMod unless they change often
  thetav <- grMod$thetav
  skel <- attr(grMod$thetav, "skel")
  thetaG <- grMod$thetaG
  thetaR <- grMod$thetaR
  nu <- grMod$nu
  sigma2e <- grMod$sigma2e
  p <- grMod$p
  lambda <- grMod$lambda
  sLc <- grMod$sLc
  AI <- grMod$AI
  dLdnu <- grMod$dLdnu


  theta <- vech2matlist(thetav, skel)
  f <- NA
  step <- 1.0
  itMat <- matrix(NA, nrow = grMod$maxit, ncol = 2 * p + 5) 
    colnames(itMat) <- c(paste0(names(thetav), "_nu"),
	paste0(names(thetav), "_theta"),
	"sigma2e", "tyPy", "logDetC", "loglik", "itTime")
  Ic <- Diagonal(x = 1, n = nrow(grMod$Cinv_ii))


  ############################################
  # 5d determine next varcomps to evaluate
  ## Evaluate and do particular REML algorithm step (EM, simplex, AI)
  #########################################################
  #########################################################
  for(i in 1:nrow(itMat)){
    vitout <- ifelse(i == 1, 0, i%%grMod$vit)
    if(grMod$v > 0 && vitout == 0){
      cat(i, "of max", grMod$maxit, "\t",
	format(Sys.time(), "%H:%M:%S"))
    }
    stItTime <- Sys.time()

    nuv <- matlist2vech(nu)
    itMat[i, 1:p] <- nuv

    if(grMod$v > 1 && vitout == 0){
      cat("\n")
      print(as.table(itMat[i, 1:p]), digits = 4, zero.print = ".")
      cat("\n")
    }



    if(lambda){
      remlOut <- reml(nu, skel, thetaG, sLc,
	grMod$modMats, grMod$W, grMod$Bpinv,
        grMod$nminffx, grMod$nminfrfx, grMod$rfxlvls, grMod$rfxIncContrib2loglik,
	thetaR = NULL,
	grMod$tWW, grMod$RHS)
    } else{
        remlOut <- reml(nu, skel, thetaG, sLc,
	  grMod$modMats, grMod$W, grMod$Bpinv,
          grMod$nminffx, grMod$nminfrfx, grMod$rfxlvls, grMod$rfxIncContrib2loglik,
	  thetaR,
	  tWW = NULL, RHS = NULL)
      }  #<-- end if/else lambda
      sigma2e[] <- remlOut$sigma2e
      grMod$sln <- remlOut$sln
      grMod$r <- remlOut$r
      sLc <- remlOut$sLc #TODO to use `update()` need to return `C` in `remlOut`

    if(lambda){
      itMat[i, (p+1):(2*p)] <- matlist2vech(nu2theta_lambda(nu, sigma2e,
							thetaG, thetaR)) 
    } else{
        itMat[i, (p+1):(2*p)] <- matlist2vech(nu2theta_noTrans(nu, thetaG,
								thetaR))
      }
    itMat[i, (2*p+1):(ncol(itMat)-1)] <- with(remlOut,
      c(sigma2e, tyPy, logDetC, loglik))

    if(grMod$v > 2 && vitout == 0){
      itMatLLcols <- match(c("sigma2e", "tyPy", "logDetC"), colnames(itMat)) 
        if(!lambda) itMatLLcols <- itMatLLcols[-1]
      print(as.table(itMat[i, itMatLLcols]), digits = 4, zero.print = ".")
      cat("\n")
    }

    # 5c check convergence criteria
    ## Knight 2008 (ch. 6) says Searle et al. 1992 and Longford 1993 discuss diff types of converg. crit.
    ## See Appendix 2 of WOMBAT help manual for 4 convergence criteria used
    cc <- rep(NA, 4)
    if(i > 1){
      # wombat 1
      cc[1] <- diff(itMat[c(i-1, i), "loglik"]) < grMod$cctol[1]
      # wombat 2 (eqn. A.1) (also Knight 2008 (eqn. 6.1) criteria
      cc[2] <- sqrt(sum((itMat[i, 1:p] - itMat[(i-1), 1:p])^2) / sum(itMat[i, 1:p]^2)) < grMod$cctol[2]
      if(grMod$algit[i] == "AI"){
        # wombat 3 (eqn. A.2): Norm of the gradient vector
        # AI only
#TODO Does this go here or maybe after AI
## Does this step happen for last AI matrix (i-1) or current (i)?
        cc[3] <- sqrt(sum(dLdnu * dLdnu)) < grMod$cctol[3]
        # wombat 4 (eqn A.3): Newton decrement (see Boyd & Vandenberghe 2004 cited in wombat)
        # AI only
#        cc[4] <- -1 * c(crossprod(dLdnu, H) %*% dLdnu)
      }
    } else cc[1] <- FALSE  #<-- ensures one of the EM/AI/etc algorithms used if i==1



    if(grMod$v > 0 && vitout == 0){
      cat("\tlL:", format(round(itMat[i, "loglik"], 6), nsmall = 6))
      if(grMod$v > 1) cat("\tConvergence crit:", cc, "\n")
    }


    if(!all(cc, na.rm = TRUE)){
      ############################
      #    EM
      ############################
      if(grMod$algit[i] == "EM"){
        if(grMod$v > 1 && vitout == 0) cat("\n\tEM to find next nu")
        emOut <- em(nuv, thetaG, thetaR,
            grMod$modMats, grMod$nminffx, sLc, grMod$ndgeninv, grMod$sln, grMod$r)
          nuvout <- emOut$nuv
          grMod$Cinv_ii <- emOut$Cinv_ii
      }


      ############################
      #    AI
      ############################
      if(grMod$algit[i] == "AI"){
        if(grMod$v > 1 && vitout == 0) cat("\n\tAI to find next nu")
#FIXME Currently, only allow when not: 
if(nrow(theta[[thetaR]]) != 1){
  stop(cat("\nAI algorithm currently only works for a single residual variance"))
}
        Cinv <- solve(a = sLc, b = Ic, system = "A")
        grMod$Cinv_ii <- diag(Cinv)

        if(lambda){
          AI <- ai(nuv, skel, thetaG,
	              grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
	  	      thetaR = NULL,
		      sigma2e)  #<-- NULL if lambda==FALSE
	  dLdnu <- gradFun(nuv, thetaG, grMod$modMats, Cinv, grMod$sln,
	    	      sigma2e = sigma2e, r = NULL, nminfrfx = NULL)
#          dLdnu_TEST <- gradFun_TEST(nuv, thetaG,
#	  	      grMod$modMats, sLc, grMod$ndgeninv, grMod$sln,	
#		      sigma2e = sigma2e,   #<-- NULL if lambda==FALSE
#		      thetaR = NULL, r = NULL, nminfrfx = NULL)  #<-- NULL if lambda==TRUE
#          dLdnu_TEST2 <- gradFun_TEST2(nuv, thetaG,
#	  	      grMod$modMats, sLc, grMod$ndgeninv, grMod$sln,	
#		      sigma2e = sigma2e,   #<-- NULL if lambda==FALSE
#		      thetaR = NULL, r = NULL, nminfrfx = NULL)  #<-- NULL if lambda==TRUE
       } else{
            AI <- ai(nuv, skel, thetaG,
        		grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
                        thetaR,   #<-- NULL if lambda==TRUE
		        sigma2e = NULL)
            dLdnu <- gradFun(nuv, thetaG, grMod$modMats, Cinv, grMod$sln,
  	      sigma2e = NULL, grMod$r, grMod$nminfrfx)
          }

        ## Find next set of parameters using a quasi-Newton method/algorithm
        ### Meyer 1989 pp. 326-327 describes quasi-Newton methods 
#TODO see Meyer 1997 eqn 58 for Marquardt 1963: theta_t+1=theta_t - (H_t + k_t * I)^{-1} g_t 
## What I do below is similar: except k_t=f
        ### Mrode 2005 eqn 11.4
        ### Johnson and Thompson 1995 eqn 12
        ####(though gremlin uses `+` instead of J & T '95 `-` because
        ##### gremlin multiplies gradient by -0.5 in `gradFun()`)

        # Check if AI can be inverted
        rcondAI <- rcond(AI)
        if(rcondAI < grMod$ezero){
          if(grMod$v > 2){
            cat("\nReciprocal condition number of AI matrix is",
	      signif(rcondAI , 2), 
	      "\n\tAI matrix may be singular - modifying diagonals")
          }  #<-- end `if v>2`
        }  #<-- end if AI singular
        ### Check/modify AI matrix to 'ensure' positive definiteness
        ### `fI` is factor to adjust AI matrix
        #### (e.g., Meyer 1997 eqn 58 and WOMBAT manual A.5 strategy 3b)
        AIeigvals <- eigen(AI, symmetric = TRUE, only.values = TRUE)$values
          d <- (3*10^-6) * AIeigvals[1]
          f <- max(0, d - AIeigvals[nrow(AI)])
        fI <- f * diag(x = 1, nrow = nrow(AI))
	##### modified 'Hessian'
        H <- fI + AI
        # Check if H can be inverted
        rcondH <- rcond(H)
        ## if H cannot be inverted do EM
        if(rcondH < grMod$ezero){
          if(grMod$v > 1){
            cat("\nReciprocal condition number of modified Hessian is",
	      signif(rcondH , 2), 
	      "\n\tHessian may be singular - switching to the EM algorithm")
          }  #<-- end `if v>1`
          if(grMod$v > 1 && vitout == 0) cat("\n\tEM to find next nu")
            emOut <- em(nuv, thetaG, thetaR,
              grMod$modMats, grMod$nminffx, sLc, grMod$ndgeninv, grMod$sln, grMod$r)
            nuvout <- emOut$nuv
            grMod$Cinv_ii <- emOut$Cinv_ii
            grMod$algit[i] <- "EM"

        } else{  #<-- end if Hessian cannot be inverted
            Hinv <- solve(H)
#TODO need a check that not proposing negative/0 variance or |correlation|>1
## Require restraining naughty components
            dnu <- Hinv %*% dLdnu   #<-- proposed change in nu parameters
 	    # Rule: if `dnu` proposed greater than 80% change in any parameter 
            ## Then implement step reduction (`grMod$step` default) else do not
            if(any(abs(dnu / matrix(nuv, ncol = 1)) > 0.8)){
	      step <- grMod$step
	    } else step <- 1.0
            nuvout <- matrix(nuv, ncol = 1) + step * dnu
            zeroV <- which(nuvout < grMod$ezero) #FIXME check variances & cov/corr separately
            if(length(zeroV) > 0L){
              if(grMod$v > 1) cat("\nVariance component(s)", zeroV, "fixed to zero")
              nuvout[zeroV] <- grMod$ezero #FIXME TODO!!!??
            }
          }  #<-- end else AI can be inverted
      }  #<-- end if algorithm is "AI"

      if(grMod$algit[i] == "bobyqa"){
stop(cat("\nNot allowing `minqa::bobyqa()` right now"))
#        if(v > 1 && vitout == 0) cat("Switching to `minqa::bobyqa()`\n")
#FIXME lower bounds if not transformed!
#        bobyout <- bobyqa(par = nuv, fn = function(x) -1*reml(x, skel), lower = ezero,
#		control = list(iprint = v, maxfun = maxit))
#        with(bobyout, cat("\t", msg, "after", feval, "iterations and with code:", ierr, "(0 desired)\n"))
#        nuout <- vech2matlist(bobyout$par, skel)
#       loglik <- -1*bobyout$fval
#FIXME do a better check of loglik and parameter changes
#	cc <- diff(c(itMat[(i-1), "loglik"], loglik)) < cctol[1] #if(bobyout$ierr == 0) TRUE else FALSE
      }
       

#TODO need to transform in order to use non-EM
##think requires obtaining gradient and hessian for both `nu` and `theta`
## See Meyer 1996 eqns ~ 45-55ish
      if(grMod$algit[i] == "NR"){
stop(cat("\nNot allowing `NR` right now"))
#        if(grMod$v > 1 && vitout == 0) cat("\n\tNR to find next nu")
#        gr <- gradFun(nuv, thetaG, thetaR, modMats, Cinv, nminfrfx, sln, r)
#        H <- hessian(func = reml, x = nuv, skel = skel) 
#tmp <- numDeriv::genD(func = reml, x = nuv, skel = skel)
#FIXME instead of `solve(H)` can I solve linear equations to give product of inverse and grad?
#        nuvout <- nuv - solve(H) %*% gr
#TODO change itnmax to correspond with algit
#tmp <- optimx(par = nuv, fn = function(x) reml(x, skel), grad = gradFun, hess = NULL,
#	lower = 0,
#	method = "L-BFGS-B",
#	itnmax = maxit, hessian = TRUE,
#	control = list(maximize = FALSE, trace = v))
      }
#        nuvout <- optim(par = nuv, fn = reml, hessian = TRUE, method = "BFGS", skel = skel)



    }  #<-- end if REML did not converge and other convergence checks







    ############################################################################
    if(lambda) nuvout[thetaR] <- nuv[thetaR]  #<-- keep R=1 (R factored out)
    nuout <- vech2matlist(nuvout, skel) 
    nu <- sapply(nuout, FUN = stTrans)
    itTime <- Sys.time() - stItTime
    if(grMod$v > 0 && vitout == 0){
      cat("\ttook ", round(itTime, 2), units(itTime), "\n")
      if(grMod$v > 2){#algit[i] == "AI" && 
        sgd <- matrix(NA, nrow = p, ncol = p+4)  #<-- `sgd` is summary.gremlinDeriv 
          dimnames(sgd) <- list(row.names(dLdnu),
            c("gradient", "", "AI", "", "AI-inv", rep("", p-1)))
        sgd[, 1] <- dLdnu
        if(grMod$algit[i] == "EM") Hinv <- H <- AI
        for(rc in 1:p){
          sgd[rc, 3:(rc+2)] <- H[rc, 1:rc]
          sgd[rc, (4+rc):(4+p)] <- Hinv[rc, rc:p]   
        }
        cat("\tstep", step, "\n")
        cat("\tAI modification", f, "\n")
        print(as.table(sgd), digits = 3, na.print = " | ", zero.print = ".")
        cat("\n")
      }  
    }
    units(itTime) <- "secs"
    itMat[i, ncol(itMat)] <- round(itTime, 1)
    if(all(cc, na.rm = TRUE)){
      if(grMod$v > 0) cat("\n***  REML converged  ***\n\n")
      break
    }

  }  # END log-likelihood iterations
  #################################### 

  itMat <- itMat[1:i, , drop = FALSE]
    rownames(itMat) <- paste(seq(i), grMod$algit[1:i], sep = "-")
  dimnames(AI) <- list(rownames(dLdnu), rownames(dLdnu))
  if(lambda){
    theta <- nu2theta_lambda(nu, sigma2e, thetaG, thetaR)
  } else{
      theta <- nu2theta_noTrans(nu, thetaG, thetaR)
    }
  thetav <- matlist2vech(theta)


  # Calculate Cinv_ii and AI for last set of parameters
  grMod$Cinv_ii <- diag(solve(a = sLc, b = Ic, system = "A"))

  ## AI
  if(lambda){
    AI <- ai(nuv, skel, thetaG,
	     grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
	     thetaR = NULL,
	     sigma2e)  #<-- NULL if lambda==FALSE

#    dLdnu <- gradFun(nuv, thetaG, grMod$modMats, Cinv, grMod$sln,
#	    	      sigma2e = sigma2e, r = NULL, nminfrfx = NULL)
  } else{
      AI <- ai(nuv, skel, thetaG,
       		grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
                thetaR,   #<-- NULL if lambda==TRUE
	        sigma2e = NULL)
#      dLdnu <- gradFun(nuv, thetaG, grMod$modMats, Cinv, grMod$sln,
#  	      sigma2e = NULL, grMod$r, grMod$nminfrfx)
      }
  

  # place these altered values back into grMod
  grMod$thetav <- thetav
  grMod$nu <- nu
  grMod$sigma2e <- sigma2e
  grMod$sLc <- sLc
  grMod$AI <- AI
  grMod$dLdnu <- dLdnu

 return(structure(list(grMod = grMod,
		itMat = itMat),
	class = "gremlinR"))
}  #<-- end `remlIt.gremlinR()`



