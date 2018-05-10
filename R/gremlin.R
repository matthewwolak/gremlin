#' Mixed-Effects REML incoporating Generalized Inverses
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
#' @import Matrix
#' @importFrom stats var
#' @importFrom methods slot
#' @importFrom methods as
#' @references
#'   Mrode. 2005.
#'   Meyer & Smith. 1996.
#'   Gilmour et al. 1995.
#'   Thompson & Johnson. 1995.
#' @seealso \code{\link[MCMCglmm:MCMCglmm-package]{MCMCglmm}}
#' @examples
#' #TODO: simple examples of the most important functions
#' \dontrun{
#'   library(nadiv)
#'   Ainv <- makeAinv(Mrode3[-c(1:2), 1:3])$Ainv
#'   mod11 <- gremlin(WWG11 ~ sex - 1,
#'	random = ~ calf,
#'	data = Mrode11,
#'	ginverse = list(calf = Ainv),
#'	Gstart = matrix(0.2), Rstart = matrix(0.4),
#'	maxit = 10, v = 2, algit = "AI")
#' }
"_PACKAGE"













################################################################################
#' Transformation of starting parameters.
#'
#' Transform start parameters into format \code{gremlin} expects.
#'
#' @aliases stTrans
#' @param x A \code{list} of starting parameters.
#' @return A sparse \sQuote{dsCMatrix}
#' @author \email{matthewwolak@@gmail.com}
stTrans <- function(x){
  if(is.numeric(x) && !is.matrix(x)) x <- as.matrix(x)
  if(!isSymmetric(x)) stop(cat("Element", x, "must be a symmetric matrix or a number\n")) 
  x <- as(x, "symmetricMatrix")
  x@uplo <- "L"
  x <- as(x, "dsCMatrix")
 x
}




################################################################################
#' Vector to list of matrices.
#'
#' Converts a vector of (co)variance parameters to a list of covariance matrices.
#'
#' @aliases vech2matlist
#' @param vech A \code{vector} of (co)variance parameters.
#' @param skeleton An example structure to map \code{vech} onto.
#' @return A list of matrices of the same structure as \code{skeleton}.
#' @author \email{matthewwolak@@gmail.com}
vech2matlist <- function(vech, skeleton){
  newmatlist <- vector("list", length = length(skeleton))
  si <- 1
  for(s in 1:length(skeleton)){
    lss <- length(skeleton[[s]][[1]])
    newmatlist[[s]] <- sparseMatrix(i = skeleton[[s]][[1]],
	p = skeleton[[s]][[2]],
	x = vech[seq(from = si, by = 1, length.out = lss)],
	dims = skeleton[[s]][[3]], symmetric = TRUE, index1 = FALSE)
    si <- si + lss
  }
  if(!is.null(names(skeleton))) names(newmatlist) <- names(skeleton)
 newmatlist
}





 
################################################################################
#' Mixed-effect modeling functions.
#'
#' Create and fit linear mixed-effect model (Gaussian data) or checking if an
#' object is a fitted model.
#'
#' @aliases gremlin gremlinR is.gremlin mkModMats
#' @param x An object of \code{class} \sQuote{gremlin}.
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
#' @param Mout A \code{logical} value that indicates if the model should quit
#'   and return the mixed model array (or M-matrix) once it is formed.
#' @param na.action What to do with NAs.
#' @param offset Should an offset be specified.
#' @param contrasts Specify the type of contrasts for the fixed effects.
#' @param Xsparse Should sparse matrices be used for the fixed effects design
#'   matrix.
#' @param \dots Additional arguments to be passed to control the model fitting.
#'
#' @return A \code{list} of class \code{gremlin} or \code{gremlinModMats}:
#'   \describe{
#'     \item{call }{The model \code{call}.}
#'     \item{modMats }{A \code{list} of the model matrices used to construct the
#'       mixed model equations.}
#'     \item{itMat }{A \code{matrix} of details about each iteration.}
#'     \item{sln }{A two column \code{matrix} of solutions and their sampling
#'       variances from the mixed model.}
#'     \item{AI }{A \code{matrix} of values containing the Average Information
#'       matrix, or second partial derivatives of the likelihood with respect to
#'       the (co)variance components. The inverse of this matrix gives the
#'       sampling variances of the (co)variance components.}
#'     \item{dLdtheta }{A single column \code{matrix} of first derivatives of
#'       the (co)variance parameters with respect to the log-Likelihood.}
#'     \item{Cinv }{The inverse coefficient \code{matrix}, in sparse matrix
#'       format, of the mixed model equations.}
#'
#'     \item{y }{The response vector.}
#'     \item{ny }{The number of responses.}
#'     \item{ncy }{The number of columns of the original response.}
#'     \item{X }{ The fixed effects design matrix.}
#'     \item{nb }{The number of columns in X.}
#'     \item{Zr }{The residual design matrix.}
#'     \item{Zg }{A list of the design matrices for each random term.}
#'     \item{nG }{The number of parameters in the G structure.}
#'     \item{listGinv }{A list of genearlized inverse matrices.}
#'     \item{logDetG }{The log-determinants of the generalized inverse matrices
#'       - necessary to calculate the log-likelihood.}
#'   }
#'
#' @references
#' Henderson
#' Mrode. 2005.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' \dontrun{
#'   library(nadiv)
#'   Ainv <- makeAinv(Mrode3[-c(1:2), 1:3])$Ainv
#'   mod11 <- gremlin(WWG11 ~ sex - 1,
#'   	random = ~ calf,
#'   	data = Mrode11,
#'   	ginverse = list(calf = Ainv),
#'   	Gstart = matrix(0.2), Rstart = matrix(0.4),
#'   	maxit = 10, v = 2, algit = "AI")
#'
#'   is(mod11)
#'
#'   mod11R <- gremlinR(WWG11 ~ sex - 1,
#'   	random = ~ calf,
#'   	data = Mrode11,
#'   	ginverse = list(calf = Ainv),
#'   	Gstart = matrix(0.2), Rstart = matrix(0.4),
#'   	maxit = 10, v = 2, algit = "AI")
#' }
#' @export
gremlin <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL,
		maxit = 20, algit = NULL,
		vit = 10, v = 1, Mout = FALSE, ...){
#FIXME USE?		control = gremlinControl(), ...){

  stopifnot(inherits(formula, "formula"), length(formula) == 3L)
  mc <- as.list(match.call())
  if(v > 0) cat("          gremlin started:\t", format(Sys.time(), "%H:%M:%S"), "\n")
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
  if(is.null(mc$algit)) algit <- c(rep("EM", 2), rep("AI", max(0, maxit-2))) else algit <- eval(mc$algit)
  if(length(algit) == 1 && algit %in% c("EM", "AI", "bobyqa")) algit <- rep(algit, maxit)
  if(is.null(mc$ezero)) ezero <- 1e-8 else ezero <- eval(mc$ezero)

  start <- list(G = Gstart, R = Rstart)
#TODO change `R.` to `R1` that way will match G1, G2, etc. for >1 G sections
##XXX then change how find thetaGorR by grep or something like it versus strsplit on `.`
  theta <- c(G = sapply(start$G, FUN = stTrans), R. = stTrans(start$R))
  thetaGorR <- sapply(strsplit(names(theta), ".", fixed = TRUE), FUN = "[[", i = 1)

#XXX Do above TODO sooner rather than later!











#FIXME ensure grep() is best way and won't mess up when multiple Gs and Rs
  thetaG <- grep("G", thetaGorR, ignore.case = FALSE)
  thetaR <- grep("R", thetaGorR, ignore.case = FALSE)
  thetav <- sapply(theta, FUN = slot, name = "x")
  skel <- lapply(seq(length(theta)), FUN = function(i){mapply(slot, theta[i], c("i", "p", "Dim"))})
    names(skel) <- names(theta)
  lthetav <- length(thetav)
#FIXME figure out environments so all in the `with(modMats...)` environment is accessible outside of it
## How to place these objects?
#  with(modMats, {



#FIXME make below uni=TRUE if R=I sigma2e
#TODO put `uni` in `mkModMats()`
    if(modMats$ncy == 1) uni <- TRUE else stop("gremlin isn't old enough to play with multivariate models")
#TODO is `start` even necessary any more? Can I just keep `Gstart` and `Rstart`

#FIXME: Add G to ginverse elements from `mkmod()` in c++ 
## to take advantage of cs_kroneckerAupdate and cs_kroneckerIupdate
#FIXME: change G to cholesky of G with log(diagonals)
## e.g., parameterisation to ensure postive-definiteness
### when to do cholesky factorization of G?
#### is it easier to do direct product of cholesky's than direct product then cholesky?
#### If the latter then save the symbolic factorization and just do updates?






    # Rand Fx incidence matrix part of 'log(|G|)'
    #FIXME: Only works for independent random effects right now!
    rfxIncContrib2loglik <- sum(unlist(modMats$logDetG))

    # 1 Create mixed model array (M) and coefficient matrix of MME (C)
    # quadratic at bottom-right (Meyer & Smith 1996, eqn 6)
    ##1a form W by `cBind()` X and each Z_i
    if(modMats$nG < 1) W <- modMats$X else  W <- cBind(modMats$X, modMats$Zg[[1]]) 
    if(modMats$nG > 1){
      for(g in 2:modMats$nG){
        W <- cBind(W, modMats$Zg[[g]])
      }
    }
    tWW <- crossprod(W)  #TODO Add statement to include `Rinv`
    zero <- Diagonal(n = modMats$nb, x = 0)  #<-- used every iteration
    #############################################################
    # REML doesn't change with any of above
    #############################################################


    # transform starting parameters to 'nu' scale
    ## cholesky of covariance matrices, then take log of transformed diagonals
    cRinv <- solve(chol(theta[[thetaR]]))
    Ginv <- lapply(thetaG, FUN = function(x){as(crossprod(cRinv, theta[[x]]) %*% cRinv, "symmetricMatrix")}) # Meyer 1991, p.77
  
    ##1c Now make coefficient matrix of MME
    ##`sapply()` to invert G_i and multiply with ginverse element (e.g., I or Ginv)
    if(modMats$nG > 0){
      C <- as(tWW + bdiag(c(zero,
	sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGinv[[u]], solve(Ginv[[u]]))}))), "symmetricMatrix")
    } else C <- as(tWW + diag(zero), "symmetricMatrix")

    ##1d Find best order of MMA/C/pivots!
    # Graser et al. 1987 (p1363) when singular C, |C| not invariant to order
    ##of C, therefore same order (of pivots) must be used in each loglik iteration
    ## Hadfield (2010, MCMCglmm paper App B): details on chol(), ordering, and updating
    # supernodal decomposition FALSE
    ## ?Cholesky documentation demos simplicial smaller sized object/more sparse
    ### However, tested `warcolak` trait1 with VA and VD univariate model
    #### `super = FALSE` faster and smaller than when `super = TRUE`
#if(any(eigen(C)$values < 0)) cat(which(eigen(C)$values < 0), "\n") #FIXME
    sLc <- Cholesky(C, perm = TRUE, LDL = FALSE, super = FALSE)
    Pc <- as(sLc, "pMatrix")
    # TODO implement options to do different methods to obtain Cinv
    ## Also, investigate/test fastest method to obtain Cinv from `sLc`
    ###XXX see note/idea in "../myNotesInsights/invFromChol.Rmd"
    Cinv <- NULL#chol2inv(sLc)
    # original order obtained by: t(Pc) %*% Lc %*% Pc or `crossprod(Pc, Lc) %*% Pc`
    D <- crossprod(modMats$y) # <-- Same every iteration
    RHSperm <- Matrix(crossprod(Pc, crossprod(W, modMats$y)), nrow = modMats$ncy, sparse = TRUE)  # <-- Same every iteration
    tRHSD <- rBind(t(RHSperm), D)
#Without permutation of C
#if(Mout) return(as(drop0(rBind(cBind(crossprod(modMats$y), crossprod(modMats$y, W)),
#	cBind(crossprod(W, modMats$y), C))), "symmetricMatrix"))
#With permutation of C
if(Mout) return(as(drop0(rBind(cBind(D, RHSperm),
	cBind(t(RHSperm), crossprod(Pc, C) %*% Pc))), "symmetricMatrix"))

# https://software.intel.com/en-us/node/531896
    ### M = [C   t(RHSperm)]   = Lm t(Lm) = [Lc                 0  ]   [t(Lc)   Lc^{-1}t(RHSperm)] 
    #       [RHSperm    D  ]                [RHSpermt(Lc^{-1})  Ls ]   [ 0      t(Ls)            ]
    ### C = Lc t(Lc)
    ### S = D - B C^{-1} t(B) the Schur complement
    ### S = Ls t(Ls)
#Alternatively: need to keep working out below
#    Lb <- tcrossprod(RHSperm, solve(sLc))
#    S <- D - tcrossprod(Lb)
#    Ls <- chol(S)
#MMA <- rBind(cBind(tcrossprod(Lc), tcrossprod(Lc, Lb)),
#	cBind(tcrossprod(Lb, Lc), tcrossprod(Lb) + tcrossprod(Ls2)))
# FIXME   tyPy <- S@x #tail(L@x, 1)^2   # faster to do `diag(L)[ncol(L)]^2`???

#TODO what to do about this: test and see which is faster (scholM TRUE/FALSE)?
## Need sLm to get tyPy? but don't use sLm otherwise below / do use C, sLc, Cinv below
#    M <- as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
#	tRHSD), "symmetricMatrix")
    # `perm=FALSE` else block diagonal no longer `bdiag(C, D)`
    sLm <- Cholesky(as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
	    tRHSD), "symmetricMatrix"), perm = FALSE, LDL = FALSE, super = FALSE)
#FIXME `tyPy` needs to be >1? So, below the `length(sLm@x)` index and/or `n=1` agrument to `tail()` needs to be changed to match dimension of y vector
    tyPy <- numeric(1)
    logDetC <- numeric(1)
    sigma2e <- numeric(1)
    # 5b log(|R|) and log(|G|) <-- Meyer 1989 (uni) & 1991 (multivar)
    # Only have to do these once per model
#TODO put these with `mkModMats()` - need to figure out multivariate version/format
    nminffx <- modMats$ny - modMats$nb
    rfxlvls <- sapply(modMats$Zg, FUN = ncol)
    nr <- if(length(rfxlvls) == 0) 0 else sum(rfxlvls)
    nminfrfx <- nminffx - nr
    sln <- matrix(0, nrow = nrow(C), ncol = 1)
    r <- matrix(0, nrow = modMats$ny, ncol = 1)

    AI <- matrix(NA, nrow = lthetav, ncol = lthetav)
    dLdtheta <- matrix(NA, nrow = lthetav, ncol = 1, dimnames = list(names(thetav), NULL))
    if("AI" %in% algit) PorVinv <- matrix(0, nrow = modMats$ny, ncol = modMats$ny)




    
    ###############################
    # reml()
    ###############################
    reml <- function(thetav, skel){
      theta <- vech2matlist(thetav, skel)
      nu <- theta
#XXX Don't need for EM algorithm
#    transform <- FALSE
#    if(transform){
#      nu <- list(G = lapply(theta$G, FUN = function(x){log(chol(x))}),
#		R = log(chol(theta$R)))
#    } else nu <- theta
      cRinv <- solve(chol(nu[[thetaR]]))
      Ginv <- lapply(thetaG, FUN = function(x){as(crossprod(cRinv, nu[[x]]) %*% cRinv, "symmetricMatrix")}) # Meyer 1991, p.77
  

      ##1c Now make coefficient matrix of MME
      ##`sapply()` to invert G_i and multiply with ginverse element (e.g., I or Ginv)
      if(modMats$nG > 0){
        C <<- as(tWW + bdiag(c(zero,
	  sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGinv[[u]], solve(Ginv[[u]]))}))), "symmetricMatrix")
      } else C <<- as(tWW + diag(zero), "symmetricMatrix")

      sLc <<- update(sLc, C)
      #TODO see note/idea in "../myNotesInsights/invFromChol.Rmd"
      #TODO should update the permutation each iteration????
      #Cinv <<- chol2inv(sLc) # chol2inv gives Cinv in same permutation as C (not permutation of sLc)
      Cinv <<- solve(C)  #<-- XXX Faster than chol2inv(sLc) atleast for warcolak
#      M <<- as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
#	    tRHSD), "symmetricMatrix")
      sLm <<- update(sLm, as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
	    tRHSD), "symmetricMatrix"))
      tyPy[] <<- sLm@x[length(sLm@x)]^2   # Meyer & Smith 1996, eqn. 14
      # above is far faster than `tail(sLm@x, 1)^2`

      # 5 record log-like, check convergence, & determine next varcomps to evaluate  
      ##5a determine log(|C|) and y'Py
      ### Meyer & Smith 1996, eqns 12-14 (and 9)
      #### Also see Meyer & Kirkpatrick 2005 GSE. eqn. 18: if cholesky of MMA = LL'
#FIXME does sLc@x[sLc@p+1] give an NA in the vector of results for Mrode11 it=1??
      logDetC[] <<- 2 * sum(log(sLc@x[sLc@p+1][1:sLc@Dim[[1]]]))
      # alternatively from sLm: `...sum(log(sLm@x[sLm@p+1][1:(sLm@Dim[[1]]-1)]))`
      ## takes longer than sLc for more complex models

      # Residual variance
      sigma2e[] <<- tyPy / nminffx

      # nminffx is tyPy/sigma2e, simplified because sigma2e = tyPy / nminffx
      loglik <- nminffx + logDetC

      # 'log(|R|)'
      #TODO: assumes X of full rank
      loglik <- loglik + nminfrfx * log(sigma2e)
      #ALTERNATIVELY: If Rinv NOT factored out of MMA `loglik <- loglik + ny * log(start$R)`

      # 'log(|G|)'
      #FIXME: Only works for independent random effects right now!
      loglik <- -0.5 * (loglik + if(modMats$nG == 0) 0 else sum(sapply(seq(modMats$nG), FUN = function(x){rfxlvls[x] * log(as.vector(Ginv[[x]]*sigma2e))})) + rfxIncContrib2loglik)
      # Below uses original starting value for residual variances - for agreement with WOMBAT
      #loglik <- -0.5 * (loglik + sum(sapply(seq(nG), FUN = function(x){rfxlvls[x] * log(as.vector(start$G[[x]]))})) + rfxIncContrib2loglik)





      # solve MME for BLUEs/BLUPs
#XXX Do I need to solve for BLUEs (see Knight 2008 for just getting BLUPs eqn 2.13-15)
      ## see Mrode 2005 chapter
      # permuted RHS of MME is `RHSperm`
      sln[] <<- solve(a = sLc, b = solve(a = sLc, b = t(RHSperm), system = "P"),
	  system = "A")
      ## Cholesky is more efficient and computationally stable
      ### see Matrix::CHMfactor-class expand note about fill-in causing many more non-zeros of very small magnitude to occur
      #### see Matrix file "CHMfactor.R" method for "determinant" (note differences with half the logdet of original matrix) and the functions:
      ##### `ldetL2up` & `destructive_Chol_update`

      # calculate residuals
      r[] <<- modMats$y - W %*% sln

     loglik

    }  #<-- end `reml()` 

#  })  #<-- end `with(modMats, ...)`









  ############################################
  # 5d determine next varcomps to evaluate
  ## Evaluate and do particular REML algorithm step (EM, simplex, AI)
  ## Meyer and Smith 1996 for algorithm using derivatives of loglik
  ### eqn 15-18 (+ eqn 33-42ish) for derivatives of tyPy and logDetC
  ### Smith 1995 for very technical details
  ## Appendix 5 of Meyer/WOMBAT help manual fro how to modify AI matrix to ensure improvements of logLike
#XXX trace of a matrix corresponds to the derivative of the determinant
  # EM refs: Hofer 1998 eqn 10-12
  ## XXX note Hofer eqn 12 missing sigma2e in last term of non-residual formula
  ### see instead Mrode 2005 (p. 241-245)
  em <- function(thetain){
    si <- modMats$nb+1
    for(g in thetaG){
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
#XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
##XXX eqn. 2.44 is the score/gradient! for a varcomp
      # current indexing faster than creating logical index matrix to multiply through Cinv for the purposes of subsetting/indexing
      thetain[[g]] <- (crossprod(sln[si:ei, , drop = FALSE], modMats$listGinv[[g]]) %*% sln[si:ei, , drop = FALSE] + tr(modMats$listGinv[[g]] %*% Cinv[si:ei, si:ei])*tail(thetav, 1)) / qi
      si <- ei+1
    }
    #FIXME make sure `nminffx` == `ncol(X)` even when reduced rank
    thetain[[thetaR]] <- crossprod(modMats$y, r) / nminffx
   thetain
  }




  ############################################
  ai <- function(thetavin){
    # Deprecated setting up `Vinv` indirectly
    #V <- kronecker(modMats$Zr, thetain[[thetaR]])
    #for(u in 1:modMats$nG){
    #  V[] <- V + with(modMats, tcrossprod(Zg[[u]] %*% kronecker(solve(listGinv[[u]]), thetain[[u]]), Zg[[u]]))
    #}
    #Vinv <- solve(V)

    lthetav <- length(thetavin)
    thetain <- vech2matlist(thetavin, skel)
    # setup Vinv directly [See Johnson & Thompson 1995, Appendix 1 (eqn A1)]
##FIXME will it work for >1 residual (co)variance????  XXX

    #TODO below removes X from W, make sure `nb` is correct when reduced rank X
    tmptmpGinv <- with(modMats, bdiag(sapply(seq(nG), FUN = function(u){kronecker(listGinv[[u]], solve(thetain[[u]]))})))
    #FIXME why the difference between using thetain versus Ginv
    #tmptmpGinv <- with(modMats, bdiag(sapply(nG, FUN = function(u){kronecker(listGinv[[u]], solve(Ginv[[u]]))})))
    Rinv <- kronecker(Diagonal(n = modMats$ny, x = 1), thetaR) #FIXME change thetaR to nu or transformed scale? 
    tZRinvZ <- with(modMats, crossprod(W[, -c(1:nb)], Rinv) %*%  W[, -c(1:nb)])
    #TODO can I use tWW below?
    PorVinv <<- with(modMats, Rinv - tcrossprod(Rinv %*% W[, -c(1:nb)] %*% solve(tZRinvZ + tmptmpGinv), W[, -c(1:nb)]) %*% Rinv)  #<-- FIXME move outside?
#FIXME why P and P2 different?
    PorVinv <<- with(modMats, PorVinv - PorVinv %*% X %*% tcrossprod(solve(crossprod(X, PorVinv) %*% X), X) %*% PorVinv)
#Cinv <- chol2inv(sLc)
#    P2 <- Rinv - tcrossprod(Rinv %*% W %*% Cinv, W) %*% Rinv #<-- See Gilmour et al. 1995 end of p.1441
#    P3 <- Diagonal(n = nrow(Rinv), x = 1) - W %*% solve(C) %*% t(W) #<-- AIreml_heritabilityPkg

    # tee = e'e
    tee <- crossprod(r)
    # trCinvGinv_gg = trace[Cinv_gg %*% Ginv_gg] | tugug = u_gg' %*% Ginv_gg %*% u_gg
    ## g is the gth component of the G-structure to model
    ## Ginv is the generalized inverse (not the inverse of the G-matrix/varcomps)
#TODO make variable `length(thetaG)`
    trCinvGinv_gg <- tugug <- as.list(rep(0, length(thetaG)))
    si <- modMats$nb+1
    for(g in thetaG){ #FIXME assumes thetaG is same length as thetavin
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      isln <- sln[si:ei, , drop = FALSE]
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      tugug[[g]] <- crossprod(isln, modMats$listGinv[[g]]) %*% isln
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      trCinvGinv_gg[[g]] <- tr(modMats$listGinv[[g]] %*% Cinv[si:ei, si:ei])
#      A <- solve(Ainv)
#      AI[g, g] <- 0.5 * (t(y) %*% P %*% A %*% P %*% A %*% P %*% y) / thetavin[g]
#FIXME Check for multivariate when theta is a matrix, but below g is assumed to be a single (co)variance
      AI[g, g] <<- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, si:ei] %*% isln) / (thetavin[g]^2))@x
#(crossprod(isln, crossprod(modMats$Zg[[g]], P)) %*% modMats$Zg[[g]] %*% isln) / (thetavin[g]^2)
      if((g+1) < lthetav){
        for(k in (g+1):(lthetav-1)){  #<-- fill upper triangle
          sk <- sum(sapply(seq(k-1), FUN = function(u){ncol(modMats$Zg[[u]])})) + modMats$nb + 1
          ek <- sk - 1 + ncol(modMats$Zg[[k]])
          AI[g, k] <<- AI[k, g] <<- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, sk:ek] %*% sln[sk:ek, , drop = FALSE]) / (thetavin[g]*thetavin[k]))@x    
        }  #<-- end 'for(k ...)`
      }  #<-- end `if()`
      AI[g, lthetav] <<- AI[lthetav, g] <<- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% r) / (thetavin[g]*thetavin[lthetav]))@x
      si <- ei+1
    }   #<-- end `for(g ...)`
    AI[lthetav, lthetav] <<- 0.5 * ((crossprod(r, PorVinv) %*% r) / thetavin[lthetav]^2)@x

    # First derivatives (gradient)
#TODO check that `nminfrfx` is still `n-p-q` when more than 1 random effect 
## ALSO check what `q` of `n-p-q` is when >1 random effect
#FIXME change `[lthetav]` below to be number of residual (co)variances
    dLdtheta[lthetav] <<- 0.5*((tee / thetavin[lthetav]^2) - (nminfrfx) / thetavin[lthetav]) 
    for(g in thetaG){
      dLdtheta[lthetav] <<- dLdtheta[lthetav] - 0.5*(trCinvGinv_gg[[g]]/thetavin[g]) 
      dLdtheta[g] <<- 0.5*(tugug[[g]]/(thetavin[g]^2) - ncol(modMats$Zg[[g]])/thetavin[g] + trCinvGinv_gg[[g]]*thetavin[lthetav]/(thetavin[g]^2))
    }

    rcondAI <- rcond(AI)
    if(rcondAI < ezero){
      if(v > 1){
        cat("Reciprocal condition number of AI matrix is", signif(rcondAI , 2), "\n\tAI matrix may be singular - switching to an iteration of the EM algorithm\n")
      }
      thetavout <- sapply(sapply(em(vech2matlist(thetavin, skel)), FUN = stTrans), FUN = slot, name = "x")
    } else{
        AIinv <- solve(AI)
#TODO need a check that not proposing negative/0 variance or |correlation|>1
## Require restraining naughty components
        thetavout <- matrix(thetavin, ncol = 1) + AIinv %*% dLdtheta
        zeroV <- which(thetavout < ezero) #FIXME check variances & cov/corr separately
        if(length(zeroV) > 0L){
          if(v > 1) cat("Variance component(s)", zeroV, "fixed to zero\n")
          thetavout[zeroV] <- ezero #FIXME TODO!!!??
        }
      }
   thetavout
  }
    




  ############################################
  gradFun <- function(thetav){
    dLdtheta <- matrix(NA, nrow = lthetav, ncol = 1)
    # tee = e'e
    tee <- crossprod(r)
    # trCinvGinv_gg = trace[Cinv_gg %*% Ginv_gg] | tugug = u_gg' %*% Ginv_gg %*% u_gg
    ## g is the gth component of the G-structure to model
    ## Ginv is the generalized inverse (not the inverse of the G-matrix/varcomps)
#TODO make variable `length(thetaG)`
    trCinvGinv_gg <- tugug <- as.list(rep(0, length(thetaG)))
    si <- modMats$nb+1
    for(g in thetaG){
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE], modMats$listGinv[[g]]) %*% sln[si:ei, drop = FALSE]
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      trCinvGinv_gg[[g]] <- tr(modMats$listGinv[[g]] %*% Cinv[si:ei, si:ei])
      si <- ei+1
    } 

    # First derivatives (gradient)
#TODO check that `nminfrfx` is still `n-p-q` when more than 1 random effect 
## ALSO check what `q` of `n-p-q` is when >1 random effect
#FIXME change `[lthetav]` below to be number of residual (co)variances
    dLdtheta[lthetav] <- 0.5*((tee / tail(thetav, 1)^2) - (nminfrfx)/tail(thetav, 1)) 
    for(g in thetaG){
      dLdtheta[lthetav] <- dLdtheta[lthetav] - 0.5*(trCinvGinv_gg[[g]]/thetav[g]) 
      dLdtheta[g] <- 0.5*(tugug[[g]]/(thetav[g]^2) - ncol(modMats$Zg[[g]])/thetav[g] + trCinvGinv_gg[[g]]*tail(thetav, 1)/(thetav[g]^2))
    }

   dLdtheta
  }





  #########################################################
  #########################################################


  itMat <- matrix(NA, nrow = maxit, ncol = lthetav+5) 
    colnames(itMat) <- c(names(thetav), "sigma2e", "tyPy", "logDetC", "loglik", "itTime")
  vitseq <- seq(0, maxit, by = vit)
  for(i in 1:nrow(itMat)){
    if(v > 0 && i %in% vitseq){
      cat("  ", i, "of max", maxit, "iterations\t",
	format(Sys.time(), "%H:%M:%S"), "\n")
    }
    stItTime <- Sys.time()
    thetav <- sapply(theta, FUN = slot, name = "x")
    loglik <- reml(thetav, skel)
    itMat[i, -ncol(itMat)] <- c(thetav, sigma2e, tyPy, logDetC, loglik) 
    # 5c check convergence criteria
    ## Knight 2008 (ch. 6) says Searle et al. 1992 and Longford 1993 discuss diff types of converg. crit.
    ## See Appendix 2 of WOMBAT help manual for 4 convergence criteria used
    cc <- rep(NA, 4)
    if(i > 1){
      # wombat 1
      cc[1] <- diff(itMat[c(i-1, i), "loglik"]) < cctol[1]
      # wombat 2 (eqn. A.1) (also Knight 2008 (eqn. 6.1) criteria
      cc[2] <- sqrt(sum((itMat[i, 1:lthetav] - itMat[(i-1), 1:lthetav])^2) / sum(itMat[i, 1:lthetav]^2)) < cctol[2]
      if(algit[i] == "AI"){
        # wombat 3 (eqn. A.2): Norm of the gradient vector
        # AI only
        cc[3] <- sqrt(sum(dLdtheta * dLdtheta)) < cctol[3]
        # wombat 4 (eqn A.3): Newton decrement (see Boyd & Vandenberghe 2004 cited in wombat)
        # AI only
#        cc[4] <- -1 * c(crossprod(dLdtheta, H) %*% dLdtheta)
      }
    }




    if(!all(cc, na.rm = TRUE)){
      if(algit[i] == "EM"){
        if(v > 1 && i %in% vitseq) cat("\t EM iteration\n")
        thetaout <- em(theta)
      }

      if(algit[i] == "AI"){
        if(v > 1 && i %in% vitseq) cat("\t AI iteration\n")
#FIXME Currently, only allow when not: 
if(nrow(theta[[thetaR]]) != 1){
  stop("AI algorithm currently only works for a single residual variance")
}
        aiout <- ai(thetav)
        thetaout <- vech2matlist(aiout, skel)
#TODO need to evaluate cc criteria now? 
## How will it work so last iteration gives AI matrix of last set of parameters and not previous
      }

      if(algit[i] == "bobyqa"){
stop("Not allowing `minqa::bobyqa()` right now")
#        if(v > 1 && i %in% vitseq) cat("Switching to `minqa::bobyqa()`\n")
#FIXME lower bounds if not transformed!
#        bobyout <- bobyqa(par = thetav, fn = function(x) -1*reml(x, skel), lower = ezero,
#		control = list(iprint = v, maxfun = maxit))
#        with(bobyout, cat("\t", msg, "after", feval, "iterations and with code:", ierr, "(0 desired)\n"))
#        thetaout <- vech2matlist(bobyout$par, skel)
#       loglik <- -1*bobyout$fval
#FIXME do a better check of loglik and parameter changes
#	cc <- diff(c(itMat[(i-1), "loglik"], loglik)) < cctol[1] #if(bobyout$ierr == 0) TRUE else FALSE
      }
       

#TODO need to transform in order to use non-EM
##think requires obtaining gradient and hessian for both `nu` and `theta`
## See Meyer 1996 eqns ~ 45-55ish
      if(algit[i] == "NR"){
        if(v > 1 && i %in% vitseq) cat("\t NR iteration\n")
#        gr <- gradFun(thetav)
#        H <- hessian(func = reml, x = thetav, skel = skel) 
#tmp <- numDeriv::genD(func = reml, x = thetav, skel = skel)
#FIXME instead of `solve(H)` can I solve linear equations to give product of inverse and grad?
#        thetavout <- thetav - solve(H) %*% gr
#        thetaout <- vech2matlist(thetavout, skel) 
#TODO change itnmax to correspond with algit
#tmp <- optimx(par = thetav, fn = function(x) reml(x, skel), grad = gradFun, hess = NULL,
#	lower = 0,
#	method = "L-BFGS-B",
#	itnmax = maxit, hessian = TRUE,
#	control = list(maximize = FALSE, trace = v))
      }
#        thetavout <- optim(par = thetav, fn = reml, hessian = TRUE, method = "BFGS", skel = skel)
    }
    theta <- sapply(thetaout, FUN = stTrans)
    itTime <- Sys.time() - stItTime
    if(v > 0 && i %in% vitseq){
      cat("\t\tlL:", format(round(loglik, 6), nsmall = 6), "\t1 iteration:",
	round(itTime, 2), units(itTime), "\n")
      if(v > 1){
#        cat("\t", colnames(itMat)[-match(c("loglik", "itTime"), colnames(itMat))], "\n", sep = "  ")
#        cat("\t", round(itMat[i, -match(c("loglik", "itTime"), colnames(itMat))], 4))
        cat("\n")
        print(as.table(itMat[i, -match(c("loglik", "itTime"), colnames(itMat))]), digits = 4, zero.print = ".")
        cat("\tConvergence crit:", cc, "\n")
      }
      if(v > 2){#algit[i] == "AI" && 
        sgd <- matrix(NA, nrow = lthetav, ncol = lthetav+4)  #<-- `sgd` is summary.gremlinDeriv 
          dimnames(sgd) <- list(row.names(dLdtheta), c("gradient", "", "AI", "", "AI-inv", rep("", lthetav-1)))
        sgd[, 1] <- dLdtheta
        AIinv <- if(algit[i] == "EM") AI else solve(AI)
        for(rc in 1:lthetav){
          sgd[rc, 3:(rc+2)] <- AI[rc, 1:rc]
          sgd[rc, (4+rc):(4+lthetav)] <- AIinv[rc, rc:lthetav]   
        }
        cat("\tAI alpha", NA, "\n") #TODO add alpha/step-halving value
        print(as.table(sgd), digits = 3, na.print = " | ", zero.print = ".")
        cat("\n")
      }  
    }
    units(itTime) <- "secs"
    itMat[i, ncol(itMat)] <- round(itTime, 1)
    if(all(cc)){
      cat("REML converged\n\n")
      break
    }

  }  # END log-likelihood iterations
  #################################### 
  itMat <- itMat[1:i, , drop = FALSE]
    rownames(itMat) <- paste(seq(i), algit[1:i], sep = "-")
  dimnames(AI) <- list(rownames(dLdtheta), rownames(dLdtheta))
#FIXME delete step: might be useful to know dimensions  if(all(is.na(AI))) AI <- NULL
#FIXME delete: useful to know dimensions, e.g., `summary.gremlin()`  if(all(is.na(dLdtheta))) dLdtheta <- NULL




#FIXME FIXME FIXME
#TODO Make sln variances come from another source (eliminate Cinv)
 return(structure(list(call = as.call(mc),
		modMats = modMats,
		itMat = itMat,
		sln = cbind(Est = c(sln), Var = diag(Cinv)),
		AI = AI, dLdtheta = dLdtheta,
		Cinv = Cinv),
	class = "gremlin"))
}
#############################
# Separating and pre-allocating P and Vinv to sparse Matrix doesn't seem to make
## much of a difference
# Also changing ginverse elements to `dsCMatrix` doesn't speedup traces, since
## they end up more or less as dense matrices but in dgCMatrix from the product
#############################
















################################################################################
#' @describeIn gremlin Version of gremlin written in R
#' @export
gremlinR <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, ginverse = NULL,
		Gstart = NULL, Rstart = NULL,
		maxit = 20, algit = NULL,
		vit = 10, v = 1, Mout = FALSE, ...){
#FIXME USE?		control = gremlinControl(), ...){

  stopifnot(inherits(formula, "formula"), length(formula) == 3L)
  mc <- as.list(match.call())
  if(v > 0) cat("          gremlin started:\t", format(Sys.time(), "%H:%M:%S"), "\n")
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
  if(is.null(mc$algit)) algit <- c(rep("EM", 2), rep("AI", max(0, maxit-2))) else algit <- eval(mc$algit)
  if(length(algit) == 1 && algit %in% c("EM", "AI", "bobyqa")) algit <- rep(algit, maxit)
  if(is.null(mc$ezero)) ezero <- 1e-8 else ezero <- eval(mc$ezero)

  start <- list(G = Gstart, R = Rstart)
#TODO change `R.` to `R1` that way will match G1, G2, etc. for >1 G sections
##XXX then change how find thetaGorR by grep or something like it versus strsplit on `.`
  theta <- c(G = sapply(start$G, FUN = stTrans), R. = stTrans(start$R))
  thetaGorR <- sapply(strsplit(names(theta), ".", fixed = TRUE), FUN = "[[", i = 1)

#XXX Do above TODO sooner rather than later!











#FIXME ensure grep() is best way and won't mess up when multiple Gs and Rs
  thetaG <- grep("G", thetaGorR, ignore.case = FALSE)
  thetaR <- grep("R", thetaGorR, ignore.case = FALSE)
  thetav <- sapply(theta, FUN = slot, name = "x")
  skel <- lapply(seq(length(theta)), FUN = function(i){mapply(slot, theta[i], c("i", "p", "Dim"))})
    names(skel) <- names(theta)
  lthetav <- length(thetav)
#FIXME figure out environments so all in the `with(modMats...)` environment is accessible outside of it
## How to place these objects?
#  with(modMats, {



#FIXME make below uni=TRUE if R=I sigma2e
#TODO put `uni` in `mkModMats()`
    if(modMats$ncy == 1) uni <- TRUE else stop("gremlin isn't old enough to play with multivariate models")
#TODO is `start` even necessary any more? Can I just keep `Gstart` and `Rstart`

#FIXME: Add G to ginverse elements from `mkmod()` in c++ 
## to take advantage of cs_kroneckerAupdate and cs_kroneckerIupdate
#FIXME: change G to cholesky of G with log(diagonals)
## e.g., parameterisation to ensure postive-definiteness
### when to do cholesky factorization of G?
#### is it easier to do direct product of cholesky's than direct product then cholesky?
#### If the latter then save the symbolic factorization and just do updates?






    # Rand Fx incidence matrix part of 'log(|G|)'
    #FIXME: Only works for independent random effects right now!
    rfxIncContrib2loglik <- sum(unlist(modMats$logDetG))

    # 1 Create mixed model array (M) and coefficient matrix of MME (C)
    # quadratic at bottom-right (Meyer & Smith 1996, eqn 6)
    ##1a form W by `cBind()` X and each Z_i
    if(modMats$nG < 1) W <- modMats$X else  W <- cBind(modMats$X, modMats$Zg[[1]]) 
    if(modMats$nG > 1){
      for(g in 2:modMats$nG){
        W <- cBind(W, modMats$Zg[[g]])
      }
    }
    tWW <- crossprod(W)  #TODO Add statement to include `Rinv`
    zero <- Diagonal(n = modMats$nb, x = 0)  #<-- used every iteration
    #############################################################
    # REML doesn't change with any of above
    #############################################################


    # transform starting parameters to 'nu' scale
    ## cholesky of covariance matrices, then take log of transformed diagonals
    cRinv <- solve(chol(theta[[thetaR]]))
    Ginv <- lapply(thetaG, FUN = function(x){as(crossprod(cRinv, theta[[x]]) %*% cRinv, "symmetricMatrix")}) # Meyer 1991, p.77
  
    ##1c Now make coefficient matrix of MME
    ##`sapply()` to invert G_i and multiply with ginverse element (e.g., I or Ginv)
    if(modMats$nG > 0){
      C <- as(tWW + bdiag(c(zero,
	sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGinv[[u]], solve(Ginv[[u]]))}))), "symmetricMatrix")
    } else C <- as(tWW + diag(zero), "symmetricMatrix")

    ##1d Find best order of MMA/C/pivots!
    # Graser et al. 1987 (p1363) when singular C, |C| not invariant to order
    ##of C, therefore same order (of pivots) must be used in each loglik iteration
    ## Hadfield (2010, MCMCglmm paper App B): details on chol(), ordering, and updating
    # supernodal decomposition FALSE
    ## ?Cholesky documentation demos simplicial smaller sized object/more sparse
    ### However, tested `warcolak` trait1 with VA and VD univariate model
    #### `super = FALSE` faster and smaller than when `super = TRUE`
#if(any(eigen(C)$values < 0)) cat(which(eigen(C)$values < 0), "\n") #FIXME
    sLc <- Cholesky(C, perm = TRUE, LDL = FALSE, super = FALSE)
    Pc <- as(sLc, "pMatrix")
    # TODO implement options to do different methods to obtain Cinv
    ## Also, investigate/test fastest method to obtain Cinv from `sLc`
    ###XXX see note/idea in "../myNotesInsights/invFromChol.Rmd"
    Cinv <- NULL#chol2inv(sLc)
    # original order obtained by: t(Pc) %*% Lc %*% Pc or `crossprod(Pc, Lc) %*% Pc`
    D <- crossprod(modMats$y) # <-- Same every iteration
    RHSperm <- Matrix(crossprod(Pc, crossprod(W, modMats$y)), nrow = modMats$ncy, sparse = TRUE)  # <-- Same every iteration
    tRHSD <- rBind(t(RHSperm), D)
#Without permutation of C
#if(Mout) return(as(drop0(rBind(cBind(crossprod(modMats$y), crossprod(modMats$y, W)),
#	cBind(crossprod(W, modMats$y), C))), "symmetricMatrix"))
#With permutation of C
if(Mout) return(as(drop0(rBind(cBind(D, RHSperm),
	cBind(t(RHSperm), crossprod(Pc, C) %*% Pc))), "symmetricMatrix"))

# https://software.intel.com/en-us/node/531896
    ### M = [C   t(RHSperm)]   = Lm t(Lm) = [Lc                 0  ]   [t(Lc)   Lc^{-1}t(RHSperm)] 
    #       [RHSperm    D  ]                [RHSpermt(Lc^{-1})  Ls ]   [ 0      t(Ls)            ]
    ### C = Lc t(Lc)
    ### S = D - B C^{-1} t(B) the Schur complement
    ### S = Ls t(Ls)
#Alternatively: need to keep working out below
#    Lb <- tcrossprod(RHSperm, solve(sLc))
#    S <- D - tcrossprod(Lb)
#    Ls <- chol(S)
#MMA <- rBind(cBind(tcrossprod(Lc), tcrossprod(Lc, Lb)),
#	cBind(tcrossprod(Lb, Lc), tcrossprod(Lb) + tcrossprod(Ls2)))
# FIXME   tyPy <- S@x #tail(L@x, 1)^2   # faster to do `diag(L)[ncol(L)]^2`???

#TODO what to do about this: test and see which is faster (scholM TRUE/FALSE)?
## Need sLm to get tyPy? but don't use sLm otherwise below / do use C, sLc, Cinv below
#    M <- as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
#	tRHSD), "symmetricMatrix")
    # `perm=FALSE` else block diagonal no longer `bdiag(C, D)`
    sLm <- Cholesky(as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
	    tRHSD), "symmetricMatrix"), perm = FALSE, LDL = FALSE, super = FALSE)
#FIXME `tyPy` needs to be >1? So, below the `length(sLm@x)` index and/or `n=1` agrument to `tail()` needs to be changed to match dimension of y vector
    tyPy <- numeric(1)
    logDetC <- numeric(1)
    sigma2e <- numeric(1)
    # 5b log(|R|) and log(|G|) <-- Meyer 1989 (uni) & 1991 (multivar)
    # Only have to do these once per model
#TODO put these with `mkModMats()` - need to figure out multivariate version/format
    nminffx <- modMats$ny - modMats$nb
    rfxlvls <- sapply(modMats$Zg, FUN = ncol)
    nr <- if(length(rfxlvls) == 0) 0 else sum(rfxlvls)
    nminfrfx <- nminffx - nr
    sln <- matrix(0, nrow = nrow(C), ncol = 1)
    r <- matrix(0, nrow = modMats$ny, ncol = 1)

    AI <- matrix(NA, nrow = lthetav, ncol = lthetav)
    dLdtheta <- matrix(NA, nrow = lthetav, ncol = 1, dimnames = list(names(thetav), NULL))
    if("AI" %in% algit) PorVinv <- matrix(0, nrow = modMats$ny, ncol = modMats$ny)




    
    ###############################
    # reml()
    ###############################
    reml <- function(thetav, skel){
      theta <- vech2matlist(thetav, skel)
      nu <- theta
#XXX Don't need for EM algorithm
#    transform <- FALSE
#    if(transform){
#      nu <- list(G = lapply(theta$G, FUN = function(x){log(chol(x))}),
#		R = log(chol(theta$R)))
#    } else nu <- theta
      cRinv <- solve(chol(nu[[thetaR]]))
      Ginv <- lapply(thetaG, FUN = function(x){as(crossprod(cRinv, nu[[x]]) %*% cRinv, "symmetricMatrix")}) # Meyer 1991, p.77
  

      ##1c Now make coefficient matrix of MME
      ##`sapply()` to invert G_i and multiply with ginverse element (e.g., I or Ginv)
      if(modMats$nG > 0){
        C <<- as(tWW + bdiag(c(zero,
	  sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGinv[[u]], solve(Ginv[[u]]))}))), "symmetricMatrix")
      } else C <<- as(tWW + diag(zero), "symmetricMatrix")

      sLc <<- update(sLc, C)
      #TODO see note/idea in "../myNotesInsights/invFromChol.Rmd"
      #TODO should update the permutation each iteration????
      #Cinv <<- chol2inv(sLc) # chol2inv gives Cinv in same permutation as C (not permutation of sLc)
      Cinv <<- solve(C)  #<-- XXX Faster than chol2inv(sLc) atleast for warcolak
#      M <<- as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
#	    tRHSD), "symmetricMatrix")
      sLm <<- update(sLm, as(cBind(rBind(crossprod(Pc, C) %*% Pc, RHSperm),
	    tRHSD), "symmetricMatrix"))
      tyPy[] <<- sLm@x[length(sLm@x)]^2   # Meyer & Smith 1996, eqn. 14
      # above is far faster than `tail(sLm@x, 1)^2`

      # 5 record log-like, check convergence, & determine next varcomps to evaluate  
      ##5a determine log(|C|) and y'Py
      ### Meyer & Smith 1996, eqns 12-14 (and 9)
      #### Also see Meyer & Kirkpatrick 2005 GSE. eqn. 18: if cholesky of MMA = LL'
#FIXME does sLc@x[sLc@p+1] give an NA in the vector of results for Mrode11 it=1??
      logDetC[] <<- 2 * sum(log(sLc@x[sLc@p+1][1:sLc@Dim[[1]]]))
      # alternatively from sLm: `...sum(log(sLm@x[sLm@p+1][1:(sLm@Dim[[1]]-1)]))`
      ## takes longer than sLc for more complex models

      # Residual variance
      sigma2e[] <<- tyPy / nminffx

      # nminffx is tyPy/sigma2e, simplified because sigma2e = tyPy / nminffx
      loglik <- nminffx + logDetC

      # 'log(|R|)'
      #TODO: assumes X of full rank
      loglik <- loglik + nminfrfx * log(sigma2e)
      #ALTERNATIVELY: If Rinv NOT factored out of MMA `loglik <- loglik + ny * log(start$R)`

      # 'log(|G|)'
      #FIXME: Only works for independent random effects right now!
      loglik <- -0.5 * (loglik + if(modMats$nG == 0) 0 else sum(sapply(seq(modMats$nG), FUN = function(x){rfxlvls[x] * log(as.vector(Ginv[[x]]*sigma2e))})) + rfxIncContrib2loglik)
      # Below uses original starting value for residual variances - for agreement with WOMBAT
      #loglik <- -0.5 * (loglik + sum(sapply(seq(nG), FUN = function(x){rfxlvls[x] * log(as.vector(start$G[[x]]))})) + rfxIncContrib2loglik)





      # solve MME for BLUEs/BLUPs
#XXX Do I need to solve for BLUEs (see Knight 2008 for just getting BLUPs eqn 2.13-15)
      ## see Mrode 2005 chapter
      # permuted RHS of MME is `RHSperm`
      sln[] <<- solve(a = sLc, b = solve(a = sLc, b = t(RHSperm), system = "P"),
	  system = "A")
      ## Cholesky is more efficient and computationally stable
      ### see Matrix::CHMfactor-class expand note about fill-in causing many more non-zeros of very small magnitude to occur
      #### see Matrix file "CHMfactor.R" method for "determinant" (note differences with half the logdet of original matrix) and the functions:
      ##### `ldetL2up` & `destructive_Chol_update`

      # calculate residuals
      r[] <<- modMats$y - W %*% sln

     loglik

    }  #<-- end `reml()` 

#  })  #<-- end `with(modMats, ...)`









  ############################################
  # 5d determine next varcomps to evaluate
  ## Evaluate and do particular REML algorithm step (EM, simplex, AI)
  ## Meyer and Smith 1996 for algorithm using derivatives of loglik
  ### eqn 15-18 (+ eqn 33-42ish) for derivatives of tyPy and logDetC
  ### Smith 1995 for very technical details
  ## Appendix 5 of Meyer/WOMBAT help manual fro how to modify AI matrix to ensure improvements of logLike
#XXX trace of a matrix corresponds to the derivative of the determinant
  # EM refs: Hofer 1998 eqn 10-12
  ## XXX note Hofer eqn 12 missing sigma2e in last term of non-residual formula
  ### see instead Mrode 2005 (p. 241-245)
  em <- function(thetain){
    si <- modMats$nb+1
    for(g in thetaG){
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
#XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
##XXX eqn. 2.44 is the score/gradient! for a varcomp
      # current indexing faster than creating logical index matrix to multiply through Cinv for the purposes of subsetting/indexing
      thetain[[g]] <- (crossprod(sln[si:ei, , drop = FALSE], modMats$listGinv[[g]]) %*% sln[si:ei, , drop = FALSE] + tr(modMats$listGinv[[g]] %*% Cinv[si:ei, si:ei])*tail(thetav, 1)) / qi
      si <- ei+1
    }
    #FIXME make sure `nminffx` == `ncol(X)` even when reduced rank
    thetain[[thetaR]] <- crossprod(modMats$y, r) / nminffx
   thetain
  }




  ############################################
  ai <- function(thetavin){
    # Deprecated setting up `Vinv` indirectly
    #V <- kronecker(modMats$Zr, thetain[[thetaR]])
    #for(u in 1:modMats$nG){
    #  V[] <- V + with(modMats, tcrossprod(Zg[[u]] %*% kronecker(solve(listGinv[[u]]), thetain[[u]]), Zg[[u]]))
    #}
    #Vinv <- solve(V)

    lthetav <- length(thetavin)
    thetain <- vech2matlist(thetavin, skel)
    # setup Vinv directly [See Johnson & Thompson 1995, Appendix 1 (eqn A1)]
##FIXME will it work for >1 residual (co)variance????  XXX

    #TODO below removes X from W, make sure `nb` is correct when reduced rank X
    tmptmpGinv <- with(modMats, bdiag(sapply(seq(nG), FUN = function(u){kronecker(listGinv[[u]], solve(thetain[[u]]))})))
    #FIXME why the difference between using thetain versus Ginv
    #tmptmpGinv <- with(modMats, bdiag(sapply(nG, FUN = function(u){kronecker(listGinv[[u]], solve(Ginv[[u]]))})))
    Rinv <- kronecker(Diagonal(n = modMats$ny, x = 1), thetaR) #FIXME change thetaR to nu or transformed scale? 
    tZRinvZ <- with(modMats, crossprod(W[, -c(1:nb)], Rinv) %*%  W[, -c(1:nb)])
    #TODO can I use tWW below?
    PorVinv <<- with(modMats, Rinv - tcrossprod(Rinv %*% W[, -c(1:nb)] %*% solve(tZRinvZ + tmptmpGinv), W[, -c(1:nb)]) %*% Rinv)  #<-- FIXME move outside?
#FIXME why P and P2 different?
    PorVinv <<- with(modMats, PorVinv - PorVinv %*% X %*% tcrossprod(solve(crossprod(X, PorVinv) %*% X), X) %*% PorVinv)
#Cinv <- chol2inv(sLc)
#    P2 <- Rinv - tcrossprod(Rinv %*% W %*% Cinv, W) %*% Rinv #<-- See Gilmour et al. 1995 end of p.1441
#    P3 <- Diagonal(n = nrow(Rinv), x = 1) - W %*% solve(C) %*% t(W) #<-- AIreml_heritabilityPkg

    # tee = e'e
    tee <- crossprod(r)
    # trCinvGinv_gg = trace[Cinv_gg %*% Ginv_gg] | tugug = u_gg' %*% Ginv_gg %*% u_gg
    ## g is the gth component of the G-structure to model
    ## Ginv is the generalized inverse (not the inverse of the G-matrix/varcomps)
#TODO make variable `length(thetaG)`
    trCinvGinv_gg <- tugug <- as.list(rep(0, length(thetaG)))
    si <- modMats$nb+1
    for(g in thetaG){ #FIXME assumes thetaG is same length as thetavin
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      isln <- sln[si:ei, , drop = FALSE]
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      tugug[[g]] <- crossprod(isln, modMats$listGinv[[g]]) %*% isln
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      trCinvGinv_gg[[g]] <- tr(modMats$listGinv[[g]] %*% Cinv[si:ei, si:ei])
#      A <- solve(Ainv)
#      AI[g, g] <- 0.5 * (t(y) %*% P %*% A %*% P %*% A %*% P %*% y) / thetavin[g]
#FIXME Check for multivariate when theta is a matrix, but below g is assumed to be a single (co)variance
      AI[g, g] <<- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, si:ei] %*% isln) / (thetavin[g]^2))@x
#(crossprod(isln, crossprod(modMats$Zg[[g]], P)) %*% modMats$Zg[[g]] %*% isln) / (thetavin[g]^2)
      if((g+1) < lthetav){
        for(k in (g+1):(lthetav-1)){  #<-- fill upper triangle
          sk <- sum(sapply(seq(k-1), FUN = function(u){ncol(modMats$Zg[[u]])})) + modMats$nb + 1
          ek <- sk - 1 + ncol(modMats$Zg[[k]])
          AI[g, k] <<- AI[k, g] <<- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, sk:ek] %*% sln[sk:ek, , drop = FALSE]) / (thetavin[g]*thetavin[k]))@x    
        }  #<-- end 'for(k ...)`
      }  #<-- end `if()`
      AI[g, lthetav] <<- AI[lthetav, g] <<- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% r) / (thetavin[g]*thetavin[lthetav]))@x
      si <- ei+1
    }   #<-- end `for(g ...)`
    AI[lthetav, lthetav] <<- 0.5 * ((crossprod(r, PorVinv) %*% r) / thetavin[lthetav]^2)@x

    # First derivatives (gradient)
#TODO check that `nminfrfx` is still `n-p-q` when more than 1 random effect 
## ALSO check what `q` of `n-p-q` is when >1 random effect
#FIXME change `[lthetav]` below to be number of residual (co)variances
    dLdtheta[lthetav] <<- 0.5*((tee / thetavin[lthetav]^2) - (nminfrfx) / thetavin[lthetav]) 
    for(g in thetaG){
      dLdtheta[lthetav] <<- dLdtheta[lthetav] - 0.5*(trCinvGinv_gg[[g]]/thetavin[g]) 
      dLdtheta[g] <<- 0.5*(tugug[[g]]/(thetavin[g]^2) - ncol(modMats$Zg[[g]])/thetavin[g] + trCinvGinv_gg[[g]]*thetavin[lthetav]/(thetavin[g]^2))
    }

    rcondAI <- rcond(AI)
    if(rcondAI < ezero){
      if(v > 1){
        cat("Reciprocal condition number of AI matrix is", signif(rcondAI , 2), "\n\tAI matrix may be singular - switching to an iteration of the EM algorithm\n")
      }
      thetavout <- sapply(sapply(em(vech2matlist(thetavin, skel)), FUN = stTrans), FUN = slot, name = "x")
    } else{
        AIinv <- solve(AI)
#TODO need a check that not proposing negative/0 variance or |correlation|>1
## Require restraining naughty components
        thetavout <- matrix(thetavin, ncol = 1) + AIinv %*% dLdtheta
        zeroV <- which(thetavout < ezero) #FIXME check variances & cov/corr separately
        if(length(zeroV) > 0L){
          if(v > 1) cat("Variance component(s)", zeroV, "fixed to zero\n")
          thetavout[zeroV] <- ezero #FIXME TODO!!!??
        }
      }
   thetavout
  }
    




  ############################################
  gradFun <- function(thetav){
    dLdtheta <- matrix(NA, nrow = lthetav, ncol = 1)
    # tee = e'e
    tee <- crossprod(r)
    # trCinvGinv_gg = trace[Cinv_gg %*% Ginv_gg] | tugug = u_gg' %*% Ginv_gg %*% u_gg
    ## g is the gth component of the G-structure to model
    ## Ginv is the generalized inverse (not the inverse of the G-matrix/varcomps)
#TODO make variable `length(thetaG)`
    trCinvGinv_gg <- tugug <- as.list(rep(0, length(thetaG)))
    si <- modMats$nb+1
    for(g in thetaG){
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE], modMats$listGinv[[g]]) %*% sln[si:ei, drop = FALSE]
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      trCinvGinv_gg[[g]] <- tr(modMats$listGinv[[g]] %*% Cinv[si:ei, si:ei])
      si <- ei+1
    } 

    # First derivatives (gradient)
#TODO check that `nminfrfx` is still `n-p-q` when more than 1 random effect 
## ALSO check what `q` of `n-p-q` is when >1 random effect
#FIXME change `[lthetav]` below to be number of residual (co)variances
    dLdtheta[lthetav] <- 0.5*((tee / tail(thetav, 1)^2) - (nminfrfx)/tail(thetav, 1)) 
    for(g in thetaG){
      dLdtheta[lthetav] <- dLdtheta[lthetav] - 0.5*(trCinvGinv_gg[[g]]/thetav[g]) 
      dLdtheta[g] <- 0.5*(tugug[[g]]/(thetav[g]^2) - ncol(modMats$Zg[[g]])/thetav[g] + trCinvGinv_gg[[g]]*tail(thetav, 1)/(thetav[g]^2))
    }

   dLdtheta
  }





  #########################################################
  #########################################################


  itMat <- matrix(NA, nrow = maxit, ncol = lthetav+5) 
    colnames(itMat) <- c(names(thetav), "sigma2e", "tyPy", "logDetC", "loglik", "itTime")
  vitseq <- seq(0, maxit, by = vit)
  for(i in 1:nrow(itMat)){
    if(v > 0 && i %in% vitseq){
      cat("  ", i, "of max", maxit, "iterations\t",
	format(Sys.time(), "%H:%M:%S"), "\n")
    }
    stItTime <- Sys.time()
    thetav <- sapply(theta, FUN = slot, name = "x")
    loglik <- reml(thetav, skel)
    itMat[i, -ncol(itMat)] <- c(thetav, sigma2e, tyPy, logDetC, loglik) 
    # 5c check convergence criteria
    ## Knight 2008 (ch. 6) says Searle et al. 1992 and Longford 1993 discuss diff types of converg. crit.
    ## See Appendix 2 of WOMBAT help manual for 4 convergence criteria used
    cc <- rep(NA, 4)
    if(i > 1){
      # wombat 1
      cc[1] <- diff(itMat[c(i-1, i), "loglik"]) < cctol[1]
      # wombat 2 (eqn. A.1) (also Knight 2008 (eqn. 6.1) criteria
      cc[2] <- sqrt(sum((itMat[i, 1:lthetav] - itMat[(i-1), 1:lthetav])^2) / sum(itMat[i, 1:lthetav]^2)) < cctol[2]
      if(algit[i] == "AI"){
        # wombat 3 (eqn. A.2): Norm of the gradient vector
        # AI only
        cc[3] <- sqrt(sum(dLdtheta * dLdtheta)) < cctol[3]
        # wombat 4 (eqn A.3): Newton decrement (see Boyd & Vandenberghe 2004 cited in wombat)
        # AI only
#        cc[4] <- -1 * c(crossprod(dLdtheta, H) %*% dLdtheta)
      }
    }




    if(!all(cc, na.rm = TRUE)){
      if(algit[i] == "EM"){
        if(v > 1 && i %in% vitseq) cat("\t EM iteration\n")
        thetaout <- em(theta)
      }

      if(algit[i] == "AI"){
        if(v > 1 && i %in% vitseq) cat("\t AI iteration\n")
#FIXME Currently, only allow when not: 
if(nrow(theta[[thetaR]]) != 1){
  stop("AI algorithm currently only works for a single residual variance")
}
        aiout <- ai(thetav)
        thetaout <- vech2matlist(aiout, skel)
#TODO need to evaluate cc criteria now? 
## How will it work so last iteration gives AI matrix of last set of parameters and not previous
      }

      if(algit[i] == "bobyqa"){
stop("Not allowing `minqa::bobyqa()` right now")
#        if(v > 1 && i %in% vitseq) cat("Switching to `minqa::bobyqa()`\n")
#FIXME lower bounds if not transformed!
#        bobyout <- bobyqa(par = thetav, fn = function(x) -1*reml(x, skel), lower = ezero,
#		control = list(iprint = v, maxfun = maxit))
#        with(bobyout, cat("\t", msg, "after", feval, "iterations and with code:", ierr, "(0 desired)\n"))
#        thetaout <- vech2matlist(bobyout$par, skel)
#       loglik <- -1*bobyout$fval
#FIXME do a better check of loglik and parameter changes
#	cc <- diff(c(itMat[(i-1), "loglik"], loglik)) < cctol[1] #if(bobyout$ierr == 0) TRUE else FALSE
      }
       

#TODO need to transform in order to use non-EM
##think requires obtaining gradient and hessian for both `nu` and `theta`
## See Meyer 1996 eqns ~ 45-55ish
      if(algit[i] == "NR"){
        if(v > 1 && i %in% vitseq) cat("\t NR iteration\n")
#        gr <- gradFun(thetav)
#        H <- hessian(func = reml, x = thetav, skel = skel) 
#tmp <- numDeriv::genD(func = reml, x = thetav, skel = skel)
#FIXME instead of `solve(H)` can I solve linear equations to give product of inverse and grad?
#        thetavout <- thetav - solve(H) %*% gr
#        thetaout <- vech2matlist(thetavout, skel) 
#TODO change itnmax to correspond with algit
#tmp <- optimx(par = thetav, fn = function(x) reml(x, skel), grad = gradFun, hess = NULL,
#	lower = 0,
#	method = "L-BFGS-B",
#	itnmax = maxit, hessian = TRUE,
#	control = list(maximize = FALSE, trace = v))
      }
#        thetavout <- optim(par = thetav, fn = reml, hessian = TRUE, method = "BFGS", skel = skel)
    }
    theta <- sapply(thetaout, FUN = stTrans)
    itTime <- Sys.time() - stItTime
    if(v > 0 && i %in% vitseq){
      cat("\t\tlL:", format(round(loglik, 6), nsmall = 6), "\t1 iteration:",
	round(itTime, 2), units(itTime), "\n")
      if(v > 1){
#        cat("\t", colnames(itMat)[-match(c("loglik", "itTime"), colnames(itMat))], "\n", sep = "  ")
#        cat("\t", round(itMat[i, -match(c("loglik", "itTime"), colnames(itMat))], 4))
        cat("\n")
        print(as.table(itMat[i, -match(c("loglik", "itTime"), colnames(itMat))]), digits = 4, zero.print = ".")
        cat("\tConvergence crit:", cc, "\n")
      }
      if(v > 2){#algit[i] == "AI" && 
        sgd <- matrix(NA, nrow = lthetav, ncol = lthetav+4)  #<-- `sgd` is summary.gremlinDeriv 
          dimnames(sgd) <- list(row.names(dLdtheta), c("gradient", "", "AI", "", "AI-inv", rep("", lthetav-1)))
        sgd[, 1] <- dLdtheta
        AIinv <- if(algit[i] == "EM") AI else solve(AI)
        for(rc in 1:lthetav){
          sgd[rc, 3:(rc+2)] <- AI[rc, 1:rc]
          sgd[rc, (4+rc):(4+lthetav)] <- AIinv[rc, rc:lthetav]   
        }
        cat("\tAI alpha", NA, "\n") #TODO add alpha/step-halving value
        print(as.table(sgd), digits = 3, na.print = " | ", zero.print = ".")
        cat("\n")
      }  
    }
    units(itTime) <- "secs"
    itMat[i, ncol(itMat)] <- round(itTime, 1)
    if(all(cc)){
      cat("REML converged\n\n")
      break
    }

  }  # END log-likelihood iterations
  #################################### 
  itMat <- itMat[1:i, , drop = FALSE]
    rownames(itMat) <- paste(seq(i), algit[1:i], sep = "-")
  dimnames(AI) <- list(rownames(dLdtheta), rownames(dLdtheta))
#FIXME delete step: might be useful to know dimensions  if(all(is.na(AI))) AI <- NULL
#FIXME delete: useful to know dimensions, e.g., `summary.gremlin()`  if(all(is.na(dLdtheta))) dLdtheta <- NULL




#FIXME FIXME FIXME
#TODO Make sln variances come from another source (eliminate Cinv)
 return(structure(list(call = as.call(mc),
		modMats = modMats,
		itMat = itMat,
		sln = cbind(Est = c(sln), Var = diag(Cinv)),
		AI = AI, dLdtheta = dLdtheta,
		Cinv = Cinv),
	class = "gremlin"))
}













################################################################################
#' @describeIn gremlinR Method to test if object is of class \code{gremlin}
#' @export
is.gremlin <- function(x) inherits(x, "gremlin")














