#' Modular mixed-effect modeling functions.
#'
#' Create and fit linear mixed-effect model (Gaussian data) in modular steps
#'
#' @aliases gremlinRmod
#' @param formula A \code{formula} for the response variable and fixed effects.
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
#' @param Gstart A \code{list} of starting (co)variance values for the the
#'   G-structure or random terms.
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
#' @param thetav,thetavin,skel variance component parameters or structure on
#' which to evaluate a model likelihood or derivatives of the likelihood with
#' respect to these parameters
#'
#' @return A \code{list} of class \code{gremlin} or \code{gremlinModMats}:
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
#'     \item{itMat }{A \code{matrix} of details about each iteration.}
#'     \item{sln }{A two column \code{matrix} of solutions and their sampling
#'       variances from the mixed model.}
#'     \item{residuals }{A \code{vector} of residual deviations, response minus
#'       the values expected based on the solutions, corresponding to the order
#'       in \code{modMats$y}.} 
#'     \item{nu }{A \code{matrix} of (co)variance components on the transformed
#'       scale at the last iteration.}
#'     \item{theta }{A \code{matrix} of (co)variance components at the last
#'       iteration.}
#'     \item{AI }{A \code{matrix} of values containing the Average Information
#'       matrix, or second partial derivatives of the likelihood with respect to
#'       the transformed (co)variance components. The inverse of this matrix
#'       gives the sampling variances of these transformed (co)variance components.}
#'     \item{dLdnu }{A single column \code{matrix} of first derivatives of
#'       the transformed (co)variance parameters with respect to the
#'       log-Likelihood.}
#'   }
#'
#' @references
#' Henderson
#' Mrode. 2005.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#'  TODO XXX XXX Make example of each modular component working with next
#'   mod11 <- gremlinRmod(WWG11 ~ sex - 1,
#'   	random = ~ calf,
#'   	data = Mrode11,
#'   	Gstart = matrix(0.2), Rstart = matrix(0.4),
#'   	maxit = 10, v = 2)
#'
#' @export
gremlinRmod <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL, Bp = NULL,
		maxit = 20, algit = NULL,
		vit = 10, v = 1, ...){
#FIXME USE?		control = gremlinControl(), ...){

  stopifnot(inherits(formula, "formula"), length(formula) == 3L)
  mc <- as.list(match.call())
  startTime <- Sys.time()
  if(v > 0) cat("gremlin started:\t\t", format(startTime, "%H:%M:%S"), "\n")
  m <- match(c("formula", "random", "rcov", "data", "subset", "ginverse", "na.action", "offset", "contrasts", "Xsparse"), names(mc), 0)
  mMmc <- as.call(c(quote(mkModMats), mc[m]))
  modMats <- eval(mMmc, parent.frame())

  #TODO check dimensions G/Rstart
#FIXME assumes univariate
  if(is.null(mc$Gstart)) Gstart <- mc$Gstart <- as.list(rep(0.1*var(modMats$y), modMats$nG))
#FIXME assumes univariate
  if(is.null(mc$Rstart)) Rstart <- mc$Rstart <- matrix(0.5*var(modMats$y))
  if(is.null(mc$maxit) && is.null(maxit)) maxit <- 20 # FIXME will this ever happen
  if(is.null(mc$vit) && is.null(vit)) vit <- 20 #FIXME will this ever happen
  # Use WOMBAT's default values for convergence criteria
  if(is.null(mc$cctol)){
    cctol <- c(5*10^-4, 10^-8, 10^-3, NULL) # [1] for AI, alternatively 10^-5 (for EM)
  } else cctol <- eval(mc$cctol)
#TODO check on validity of inputted algorithms (format and matching to actual ones)
## Otherwise, get obscure warning about not finding `thetaout`/`nuout`
## Something like the following line, but implement partial matching
  if(!all(unique(algit) %in% c("EM", "AI", "bobyqa", "NR"))){  #TODO Update choices of algorithm if add/subtract any
    stop(cat("Algorithms:", unique(algit)[which(!unique(algit) %in% c("EM", "AI", "bobyqa", "NR"))],
      "not a valid choice. Please check values given to the `algit` argument\n"))
  }
  if(is.null(mc$algit)) algit <- c(rep("EM", 2), rep("AI", max(0, maxit-2))) else algit <- eval(mc$algit)
  if(length(algit) == 1 && algit %in% c("EM", "AI", "bobyqa")) algit <- rep(algit, maxit)
  if(is.null(mc$ezero)) ezero <- 1e-8 else ezero <- eval(mc$ezero)

#TODO change `R.` to `R1` that way will match G1, G2, etc. for >1 G sections
##XXX then change how find thetaGorR by grep or something like it versus strsplit on `.`
  theta <- c(G = sapply(Gstart, FUN = stTrans), R. = stTrans(Rstart))
  thetaGorR <- sapply(strsplit(names(theta), ".", fixed = TRUE), FUN = "[[", i = 1)

#XXX Do above TODO sooner rather than later!











#FIXME ensure grep() is best way and won't mess up when multiple Gs and **Rs**
  thetaG <- grep("G", thetaGorR, ignore.case = FALSE)
  thetaR <- grep("R", thetaGorR, ignore.case = FALSE)
  thetav <- sapply(theta, FUN = slot, name = "x")
  skel <- lapply(seq(length(theta)), FUN = function(i){mapply(slot, theta[i], c("i", "p", "Dim"))})
    names(skel) <- names(theta)
  p <- length(thetav)










#TODO TODO TODO TODO TODO TODO TODO
#XXX Figure out how to identify `lambda` model from model call

#FIXME quick fix for now
lambda <- length(thetaR) == 1
#TODO quick fix to turn this off from the call if I don't want to do lambda model
#TODO make `lambda == FALSE` if EM is part of algit (not sure if EM works for lambda or not) #<-- FIXME figure out how to do EM on lambda







#FIXME figure out environments so all in the `with(modMats...)` environment is accessible outside of it
## How to place these objects?
#  with(modMats, {



#FIXME make below uni=TRUE if R=I sigma2e
#TODO put `uni` in `mkModMats()`
    if(modMats$ncy == 1) uni <- TRUE else stop("gremlin isn't old enough to play with multivariate models")




#FIXME: change G to cholesky of G with log(diagonals)
## e.g., parameterisation to ensure postive-definiteness
### when to do cholesky factorization of G?
#### is it easier to do direct product of cholesky's than direct product then cholesky?
#### If the latter then save the symbolic factorization and just do updates?
#XXX Don't need for EM algorithm
  transform <- FALSE
  if(transform){ 
    #TODO FIXME need to take log of just diagonals
    #TODO also convert to lambda scale
    nu <- list(G = lapply(theta$G, FUN = function(x){log(chol(x))}),
        	R = log(chol(theta$R)))
    #TODO FIXME when converting back, do exp() of just diagonals
    #TODO also back-convert from lambda scale
    nu2theta <- function(nu){
	list(G = lapply(nu[[1L]], FUN = function(x){exp(crossprod(x))}),
	     R = exp(crossprod(nu[[2L]])))
    }
  } else{
      nu <- theta
      if(lambda){
        cRinv <- solve(chol(theta[[thetaR]]))
        nu[thetaG] <- lapply(thetaG, FUN = function(x){as(crossprod(cRinv, nu[[x]]) %*% cRinv, "symmetricMatrix")}) # Meyer 1991, p.77 (<-- ?p70/eqn4?)
        #TODO what to do if R is a matrix
        nu[[thetaR]] <- as(crossprod(cRinv, nu[[thetaR]]) %*% cRinv, "symmetricMatrix")
        nu2theta <- function(nu, sigma2e){ #TODO FIXME
          cR <- chol(matrix(sigma2e))
          theta <- nu
          theta[thetaG] <- lapply(thetaG, FUN = function(x){as(crossprod(cR, nu[[x]]) %*% cR, "symmetricMatrix")})
          theta[[thetaR]] <- as(crossprod(cR, nu[[thetaR]]) %*% cR, "symmetricMatrix")
         return(theta)
        }  #<-- end `if lambda`

      } else{
          nu2theta <- function(nu) nu #TODO FIXME
        }  #<-- end `if NOT lambda`
    }





    # 1 Create coefficient matrix of MME (C)
    ##1a form W by `cbind()` X and each Z_i
    if(modMats$nG < 1) W <- modMats$X else  W <- cbind(modMats$X, modMats$Zg[[1]]) 
    if(modMats$nG > 1){
      for(g in 2:modMats$nG){
        W <- cbind(W, modMats$Zg[[g]])
      }
    }
    # Rand Fx incidence matrix part of 'log(|G|)'
    #FIXME: Only works for independent random effects right now!
    rfxIncContrib2loglik <- sum(unlist(modMats$logDetG))
    if(is.null(mc$Bp) && is.null(Bp)){
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
    # Find Non-diagonal ginverses
    ndgeninv <- sapply(seq(modMats$nG),
	FUN = function(g){class(modMats$listGeninv[[g]]) != "ddiMatrix"})  #0=I; 1=A 
    dimsZg <- sapply(seq(modMats$nG),
	FUN = function(g){slot(modMats$Zg[[g]], "Dim")})
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
    rfxlvls <- sapply(modMats$Zg, FUN = ncol)
    nr <- if(length(rfxlvls) == 0) 0 else sum(rfxlvls)
    nminfrfx <- nminffx - nr

    AI <- matrix(NA, nrow = p, ncol = p)
    f <- NA
    dLdnu <- matrix(NA, nrow = p, ncol = 1, dimnames = list(names(thetav), NULL))
    sigma2e <- if(lambda) numeric(1) else NA #<-- only needed for `lambda` model


  itMat <- matrix(NA, nrow = maxit, ncol = 2 * p + 5) 
    colnames(itMat) <- c(paste0(names(thetav), "_nu"),
	paste0(names(thetav), "_theta"),
	"sigma2e", "tyPy", "logDetC", "loglik", "itTime")
    #############################################################
    # REML doesn't change with any of above
    #############################################################




#  })  #<-- end `with(modMats, ...)`






  ############################################
  # 5d determine next varcomps to evaluate
  ## Evaluate and do particular REML algorithm step (EM, simplex, AI)
  #########################################################
  #########################################################
  for(i in 1:nrow(itMat)){
    vitout <- ifelse(i == 1, 0, i%%vit)
    if(v > 0 && vitout == 0){
      cat("  ", i, "of max", maxit, "\t\t\t",
	format(Sys.time(), "%H:%M:%S"), "\n")
    }
    stItTime <- Sys.time()

    if(lambda){
      remlOut <- reml(nu, skel, thetaG, sLc,
	modMats, W, Bpinv, nminffx, nminfrfx, rfxlvls, rfxIncContrib2loglik,
	thetaR = NULL,
	tWW, RHS)
    } else{
        remlOut <- reml(nu, skel, thetaG, sLc,
	  modMats, W, Bpinv, nminffx, nminfrfx, rfxlvls, rfxIncContrib2loglik.
	  thetaR,
	  tWW = NULL, RHS = NULL)
      }  #<-- end if/else lambda
      nuv <- sapply(nu, FUN = slot, name = "x")
      sigma2e[] <- remlOut$sigma2e
      sln <- remlOut$sln
      r <- remlOut$r
      sLc <- remlOut$sLc #TODO to use `update()` need to return `C` in `remlOut`

    itMat[i, -ncol(itMat)] <- c(nuv,
      sapply(if(lambda) nu2theta(nu, sigma2e) else nu2theta(nu),
            FUN = slot, name = "x"),
      sigma2e, remlOut$tyPy, remlOut$logDetC, remlOut$loglik) 
    # 5c check convergence criteria
    ## Knight 2008 (ch. 6) says Searle et al. 1992 and Longford 1993 discuss diff types of converg. crit.
    ## See Appendix 2 of WOMBAT help manual for 4 convergence criteria used
    cc <- rep(NA, 4)
    if(i > 1){
      # wombat 1
      cc[1] <- diff(itMat[c(i-1, i), "loglik"]) < cctol[1]
      # wombat 2 (eqn. A.1) (also Knight 2008 (eqn. 6.1) criteria
      cc[2] <- sqrt(sum((itMat[i, 1:p] - itMat[(i-1), 1:p])^2) / sum(itMat[i, 1:p]^2)) < cctol[2]
      if(algit[i] == "AI"){
        # wombat 3 (eqn. A.2): Norm of the gradient vector
        # AI only
        cc[3] <- sqrt(sum(dLdnu * dLdnu)) < cctol[3]
        # wombat 4 (eqn A.3): Newton decrement (see Boyd & Vandenberghe 2004 cited in wombat)
        # AI only
#        cc[4] <- -1 * c(crossprod(dLdnu, H) %*% dLdnu)
      }
    } else cc[1] <- FALSE  #<-- ensures one of the EM/AI/etc algorithms used if i==1





    if(!all(cc, na.rm = TRUE)){
      if(algit[i] == "EM"){
        if(v > 1 && vitout == 0) cat("\tEM to find next theta\n")
        emOut <- em(nuv, thetaG, thetaR,
		modMats, nminffx, sLc, ndgeninv, sln, r)
          nuvout <- emOut$nuv
          Cinv_ii <- emOut$Cinv_ii
        nuout <- vech2matlist(nuvout, skel) 
      }

      if(algit[i] == "AI"){
        if(v > 1 && vitout == 0) cat("\tAI to find next theta\n")
#FIXME Currently, only allow when not: 
if(nrow(theta[[thetaR]]) != 1){
  stop("AI algorithm currently only works for a single residual variance")
}
        Cinv <- solve(a = sLc, b = Ic, system = "A")
        Cinv_ii <- diag(Cinv)
        if(lambda){
          aiout <- ai(nuv, skel, thetaG,
	              modMats, W, sLc, sln, r,
	  	      thetaR = NULL,
		      sigma2e)  #<-- non-NULL if lambda==TRUE
          dLdnu <- gradFun(nuv, thetaG, modMats, Cinv, sln,
	    sigma2e = sigma2e, r = NULL, nminfrfx = NULL)
        } else{
            aiout <- ai(nuv, skel, thetaG,
        		modMats, W, sLc, sln, r,
	  		     thetaR,   #<-- non-NULL if lambda==FALSE
		             sigma2e = NULL)
            dLdnu <- gradFun(nuv, thetaG, modMats, Cinv, sln,
  	      sigma2e = NULL, r, nminfrfx)
          }
        AI <- aiout$AI
	AIinv <- solve(AI)

        ## Find next set of parameters using a quasi-Newton method/algorithm
        ### Meyer 1989 pp. 326-327 describes quasi-Newton methods 
#TODO see Meyer 1997 eqn 58 for Marquardt 1963: theta_t+1=theta_t - (H_t + k_t * I)^{-1} g_t 
        ### Mrode 2005 eqn 11.4
        ### Johnson and Thompson 1995 eqn 12
        ####(though gremlin uses `+` instead of J & T '95 `-` because
        ##### gremlin multiplies gradient by -0.5 in `gradFun()`)

        ### Check/modify AI matrix to 'ensure' positive definiteness
        ### `fI` is factor to adjust AI matrix
        #### (e.g., Meyer 1997 eqn 58 and WOMBAT manual A.5 strategy 3b)
        AIeigvals <- eigen(AI)$values
          d <- (3*10^-6) * AIeigvals[1]
          f <- max(0, d - AIeigvals[nrow(AI)])
        fI <- f * diag(x = 1, nrow = nrow(AI))
	##### modified 'Hessian'
        H <- fI + AI
        # Check if AI can be inverted
        rcondH <- rcond(H)
        ## if AI cannot be inverted do EM
        if(rcondH < ezero){
          if(v > 1){
            cat("Reciprocal condition number of AI matrix is", signif(rcondH , 2), "\n\tAI matrix may be singular - switching to an iteration of the EM algorithm\n")
          }  #<-- end `if v>1`
          if(v > 1 && vitout == 0) cat("\tEM to find next theta\n")
            emOut <- em(nuv, thetaG, thetaR,
		modMats, nminffx, sLc, ndgeninv, sln, r)
            nuvout <- emOut$nuv
            Cinv_ii <- emOut$Cinv_ii
        } else{  #<-- end if AI cannot be inverted
            Hinv <- solve(H)
#TODO need a check that not proposing negative/0 variance or |correlation|>1
## Require restraining naughty components
            nuvout <- matrix(nuv, ncol = 1) + Hinv %*% dLdnu
            zeroV <- which(nuvout < ezero) #FIXME check variances & cov/corr separately
            if(length(zeroV) > 0L){
              if(v > 1) cat("Variance component(s)", zeroV, "fixed to zero\n")
              nuvout[zeroV] <- ezero #FIXME TODO!!!??
            }
          }  #<-- end else AI can be inverted

        if(lambda) nuvout[thetaR] <- nuv[thetaR]  #<-- keep R=1 (R factored out)
        nuout <- vech2matlist(nuvout, skel) 
      }  #<-- end if algorithm is "AI"

      if(algit[i] == "bobyqa"){
stop("Not allowing `minqa::bobyqa()` right now")
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
      if(algit[i] == "NR"){
        if(v > 1 && vitout == 0) cat("\tNR to find next theta\n")
#        gr <- gradFun(nuv, thetaG, thetaR, modMats, Cinv, nminfrfx, sln, r)
#        H <- hessian(func = reml, x = nuv, skel = skel) 
#tmp <- numDeriv::genD(func = reml, x = nuv, skel = skel)
#FIXME instead of `solve(H)` can I solve linear equations to give product of inverse and grad?
#        nuvout <- nuv - solve(H) %*% gr
#        nuout <- vech2matlist(nuvout, skel) 
#TODO change itnmax to correspond with algit
#tmp <- optimx(par = nuv, fn = function(x) reml(x, skel), grad = gradFun, hess = NULL,
#	lower = 0,
#	method = "L-BFGS-B",
#	itnmax = maxit, hessian = TRUE,
#	control = list(maximize = FALSE, trace = v))
      }
#        nuvout <- optim(par = nuv, fn = reml, hessian = TRUE, method = "BFGS", skel = skel)
#        nuout <- vech2matlist(nuvout, skel) 
    }


    ############################################################################
    nu <- sapply(nuout, FUN = stTrans)
    itTime <- Sys.time() - stItTime
    if(v > 0 && vitout == 0){
      cat("\tlL:", format(round(itMat[i, "loglik"], 6), nsmall = 6), "\ttook ",
	round(itTime, 2), units(itTime), "\n")
      if(v > 1){
        cat("\n")
        print(as.table(itMat[i, -match(c("loglik", "itTime"), colnames(itMat))]), digits = 4, zero.print = ".")
        cat("\tConvergence crit:", cc, "\n")
      }
      if(v > 2){#algit[i] == "AI" && 
        sgd <- matrix(NA, nrow = p, ncol = p+4)  #<-- `sgd` is summary.gremlinDeriv 
          dimnames(sgd) <- list(row.names(dLdnu), c("gradient", "", "AI", "", "AI-inv", rep("", p-1)))
        sgd[, 1] <- dLdnu
        if(algit[i] == "EM") AIinv <- AI
        for(rc in 1:p){
          sgd[rc, 3:(rc+2)] <- AI[rc, 1:rc]
          sgd[rc, (4+rc):(4+p)] <- AIinv[rc, rc:p]   
        }
#        cat("\tAI alpha", NA, "\n") #TODO add alpha/step-halving value
        cat("\tAI modification", f, "\n")
        print(as.table(sgd), digits = 3, na.print = " | ", zero.print = ".")
        cat("\n")
      }  
    }
    units(itTime) <- "secs"
    itMat[i, ncol(itMat)] <- round(itTime, 1)
    if(all(cc, na.rm = TRUE)){
      cat("REML converged\n\n")
      break
    }

  }  # END log-likelihood iterations
  #################################### 

  itMat <- itMat[1:i, , drop = FALSE]
    rownames(itMat) <- paste(seq(i), algit[1:i], sep = "-")
  dimnames(AI) <- list(rownames(dLdnu), rownames(dLdnu))
  theta <- if(lambda) nu2theta(nu, sigma2e) else nu2theta(nu)
    thetav <- sapply(theta, FUN = slot, name = "x") 


#TODO calculate AI/AIinv for last set of parameters (sampling covariances for final varcomps
###TODO make sure final varcomps returned are ones for which AI was evaluated (if using AI from above)

 endTime <- Sys.time()
 if(v > 0) cat("gremlin ended:\t\t", format(endTime, "%H:%M:%S"), "\n")





#TODO fix return values:
## are sln same as on theta scale?
## are Cinv_ii (var of sln) same as on theta scale
 return(structure(list(call = as.call(mc),
		modMats = modMats,
		itMat = itMat,
		sln = as(cbind(Est = sln, Var = Cinv_ii), "matrix"),
		residuals = as(r, "matrix"),
		nu = matrix(nuv, nrow = p, ncol = 1,
		  dimnames = list(paste0(names(thetav), "_nu"), NULL)),
		theta = matrix(thetav, nrow = p, ncol = 1,
		  dimnames = list(names(thetav), NULL)),
		AI = AI, dLdnu = dLdnu),
	class = "gremlin",
	startTime = startTime, endTime = endTime))
}  #<-- end `gremlinRmod()`


