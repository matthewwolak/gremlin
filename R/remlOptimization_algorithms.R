###############################
# reml()
###############################
#' REML optimization algorithms for mixed-effect models.
#'
#' Evaluate the REML likelihood and algorithms for iterating to find maximum 
#'   REML estimates.
#'
#' @aliases reml ai em gradFun
#' @param nu,nuvin A \code{list} or \code{vector} of (co)variance parameters to
#'   estimate on the transformed, or nu, scale.
#' @param skel A skeleton for reconstructing the list of (co)variance parameters.
#' @param thetaG,thetaR \code{Integer} vectors indexing the G-structure or
#'   R-structure of the (co)variance parameters.
#' @param sigma2e A \code{numeric} value for the residual variance estimate
#'   when it has been factored out of the Coefficient matrix of the Mixed Model
#'   Equations, thus converting the (co)variance components to ratios
#'   (represented by the variable lambda).
#' @param conv A \code{character} vector of (co)variance parameter constraints.
#' @param sLc A sparse \code{Matrix} containing the symbolic Cholesky
#'   factorization of the coefficient matrix of the Mixed Model Equations.
#' @param Cinv A sparse \code{Matrix} containing the inverse of the Coefficient
#'   matrix to the Mixed Model Equations.
#' @param modMats A \code{list} of the model matrices used to construct the
#'   mixed model equations.
#' @param W,tWW A sparse \code{Matrix} containing the design matrices for the fixed
#'   and random effects (\code{W}) and the cross-product of this (\code{tWW}).
#' @param Bpinv A matrix inverse of the matrix containing the prior specification
#'   for fixed effects.
#' @param nminffx,nminfrfx,rfxlvls \code{Integers} specifying: (1) the difference
#'   between the number of observations and fixed effects (of the full rank fixed
#'   effects design matrix (X), (2) \code{nminffx} minus the total number of
#'   random effects, and (3) a \code{vector} of levels for each term in the 
#'   random effects.
#' @param rfxIncContrib2loglik A \code{numeric} indicating the sum of constraint
#'   contributions to the log-likelihood across all terms in the random effects
#'   that have non-diagonal generalized inverse matrices (\code{ginverse}).
#'   associated with a generalized inverse (\code{ginverse}).
#' @param RHS A sparse \code{Matrix} containing the Right-Hand Side to the 
#'   Mixed Model Equations.
#' @param sln,r Sparse \code{Matrices} containing the solutions or residuals
#'   of the Mixed Model Equations.
#'
#' @return A \code{list} or \code{matrix} containing any of the previous
#'   parameters described above, or the following that are in addition to or
#'   instead of parameters above:
#'   \describe{
#'     \item{loglik }{The REML log-likelihood.}
#'     \item{tyPy,logDetC }{Components of the REML log-likelihood derived from the 
#'       Cholesky factor of the Coefficient matrix to the Mixed Model Equations.}
#'     \item{Cinv_ii }{A vector containing the diagonal elements of the inverse
#'       of the Coefficient matrix to the Mixed Model Equations (i.e., the
#'       diagonal entries of \code{Cinv}).}
#'     \item{AI }{A \code{matrix} of values containing the Average Information
#'       matrix, or second partial derivatives of the likelihood with respect to
#'       the transformed (co)variance components (nu). The inverse of this matrix
#'       gives the sampling variances of these transformed (co)variance components.}
#'     \item{dLdnu }{A single column \code{matrix} of first derivatives of
#'       the transformed (co)variance parameters (nu) with respect to the
#'       log-Likelihood.}
#'   }
#'
#' @author \email{matthewwolak@@gmail.com}
#' @export
reml <- function(nu, skel, thetaG, sLc,
	modMats, W, Bpinv, nminffx, nminfrfx, rfxlvls, rfxIncContrib2loglik,
	thetaR = NULL,  #<-- non-NULL if lambda==FALSE
	tWW = NULL, RHS = NULL){  #<-- non-NULL if lambda==TRUE

  lambda <- is.null(thetaR)
  #FIXME `[[thetaG+1]]` is just a kluge below: check when/if have > R matrices
  Rinv <- as(solve(nu[[length(thetaG)+1]]), "symmetricMatrix")
  Ginv <- lapply(thetaG, FUN = function(x){as(solve(nu[[x]]), "symmetricMatrix")}) # Meyer 1991, p.77 (<-- also see MMA on p70/eqn4)

  ##1c Now make coefficient matrix of MME
  if(lambda){
    tWKRinvW <- tWW  #<-- if lambda, data "quadratic" with R factored out
    tyRinvy <- crossprod(modMats$y)
  } else{
      ### Rinv Kronecker with Diagonal/I
      KRinv <- kronecker(Rinv, Diagonal(x = 1, n = modMats$Zr@Dim[[2L]]))
      tyRinvy <- as(crossprod(modMats$y, KRinv) %*% modMats$y, "sparseMatrix")
        ### Components of Meyer 1989 eqn 2
        tWKRinv <- crossprod(W, KRinv)
        tWKRinvW <- tWKRinv %*% W  
      RHS <- Matrix(tWKRinv %*% modMats$y, sparse = TRUE)
       ## Meyer '97 eqn 11
       ### (Note different order of RHS from Meyer '89 eqn 6; Meyer '91 eqn 4)
    }

  ##`sapply()` to multiply inverse G_i with ginverse element (e.g., I or geninv)
  if(modMats$nG > 0){
    C <- as(tWKRinvW + bdiag(c(Bpinv,
      sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGeninv[[u]], Ginv[[u]])}))), "symmetricMatrix")
  } else C <- as(tWKRinvW + Bpinv, "symmetricMatrix")
  ## Make/Update symbolic factorization with variance parameters from this iteration
  if(is.null(sLc)){
    ##1d Find best order of MMA/C/pivots!
    # Graser et al. 1987 (p1363) when singular C, |C| not invariant to order
    ##of C, therefore same order (of pivots) must be used in each loglik iteration
    ## Hadfield (2010, MCMCglmm paper App B): details on chol(), ordering, and updating
    # supernodal decomposition FALSE
    ## ?Cholesky documentation demos simplicial smaller sized object/more sparse
    ### However, tested `warcolak` trait1 with VA and VD univariate model
    #### `super = FALSE` faster and smaller than when `super = TRUE`
#if(any(eigen(C)$values < 0)) cat(which(eigen(C)$values < 0), "\n") #FIXME
    #TODO test if super=TRUE better
    ## supernodal might be better as coefficient matrix (C) becomes more dense
    #### e.g., with a genomic relatedness matrix (see Masuda et al. 2014 & 2015)
    sLc <- Cholesky(C, perm = TRUE, LDL = FALSE, super = FALSE)
    # original order obtained by: t(P) %*% L %*% P or `crossprod(P, L) %*% P`
  } else sLc <- update(sLc, C)



  # solve MME for BLUEs/BLUPs
  ## Do now, because needed as part of calculating the log-likelihood
#XXX Do I need to solve for BLUEs (see Knight 2008 for just getting BLUPs eqn 2.13-15)
  # chol2inv: Cinv in same permutation as C (not permutation of sLc/sLm)
  #Cinv <<- chol2inv(sLc) #<-- XXX ~10x slower than `solve(C)` atleast for warcolak
#      Cinv <<- solve(a = sLc, b = Ic, system = "A") #<-- XXX comparable speed to `solve(C)` at least for warcolak
  #Cinv <<- solve(C)
  ##XXX NOTE Above Cinv is in original order of C and NOT permutation of M

  sln <- solve(a = sLc, b = RHS, system = "A")
  ## Cholesky is more efficient and computationally stable
  ### see Matrix::CHMfactor-class expand note about fill-in causing many more non-zeroes of very small magnitude to occur
  #### see Matrix file "CHMfactor.R" method for "determinant" (note differences with half the logdet of original matrix) and the functions:
  ##### `ldetL2up` & `destructive_Chol_update`



  # 5 record log-like, check convergence, & determine next varcomps to evaluate  
  ##5a determine log(|C|) and y'Py
  # Boldman and Van Vleck 1991 eqn. 6 (tyPy) and 7 (log|C|)    
  # Meyer 1997 eqn. 13 (tyPy)
  tyPy <- tyRinvy - crossprod(sln, RHS)
  logDetC <- 2 * sum(log(sLc@x[sLc@p+1][1:sLc@Dim[[1L]]]))
  # alternatively see `determinant` method for CHMfactor

  # Factored out residual variance (only for lambda scale)
  sigma2e <- if(lambda) tyPy / nminffx else NA

  # Construct the log-likelihood (Meyer 1997, eqn. 8)
  ## (first put together as `-2log-likelihood`)
  ## 'log(|R|)'
  ### assumes X of full rank
  if(lambda){
    loglik <- nminfrfx * log(sigma2e)
  } else{
      #  Meyer 1997 eqn. 10
      loglik <- modMats$ny * log(nu[[thetaR]])
    }

  ## 'log(|G|)'
  #FIXME: Only works for independent random effects right now!
  if(lambda){
    logDetGfun <- function(x){rfxlvls[x] * log(as.vector(nu[[x]]*sigma2e))}
  } else{
      # Meyer 1997 eqn. 9
      logDetGfun <- function(x){rfxlvls[x] * log(as.vector(nu[[x]]))}
    }
  if(modMats$nG != 0){
    loglik <- loglik + sum(sapply(seq(modMats$nG), FUN = logDetGfun)) + rfxIncContrib2loglik
  }

  ## log(|C|) + tyPy
  ### if `lambda` then nminffx=tyPy/sigma2e, simplified because sigma2e=tyPy/nminffx
  ### and mulitply by -0.5 to calculate `loglik` from `-2loglik`
  loglik <- -0.5 * (loglik + logDetC + if(lambda) nminffx else tyPy)




  # calculate residuals
  r <- modMats$y - W %*% sln

 return(list(loglik = loglik@x,
		sigma2e = if(lambda) sigma2e@x else NA,
		tyPy = tyPy@x, logDetC = logDetC,
		sln = sln, r = r, sLc = sLc))
}  #<-- end `reml()` 


################################################################################



















################################################################################
## Meyer and Smith 1996 for algorithm using derivatives of loglik
### eqn 15-18 (+ eqn 33-42ish) for derivatives of tyPy and logDetC
### Smith 1995 for very technical details
# EM refs: Hofer 1998 eqn 10-12
## note Hofer eqn 12 has no sigma2e in last term of non-residual formula
### ?mistake? in Mrode 2005 (p. 241-245), which does include sigma2e
#' @rdname reml
#' @param tugug A list of numeric values for the $u_g' u_g$ products of the
#'   solution vector for each of the \code{g} variance parameters. 
#' @param trace A list of traces of the inverse coefficient matrix for each of
#'   of the \code{g} variance parameters.
#' @param y The response vector.
#' @export
em <- function(nuvin, thetaG, thetaR, conv,
	rfxlvls, tugug, trace, y = NULL, r = NULL, nminffx = NULL){

  # if no variance components (except residual), skip to the end
  if(length(thetaG) > 0){
    for(g in thetaG){
      if(conv[g] == "F") next  #<-- skip if parameter Fixed
      qi <- rfxlvls[g]

      #TODO check `*tail(nuv,1)` correctly handles models with covariance matrices
      nuvin[g] <- as(as(matrix((tugug[[g]] + trace[[g]]) / qi),  #XXX was `trace*tail(nuvin,1)` 
          "symmetricMatrix"),
          "dsCMatrix")
    }  #<-- end `for g`
  }  #<-- end if no variance components besides residuals

  if(conv[thetaR] != "F") nuvin[thetaR] <- crossprod(y, r) / nminffx

 return(nuvin)
}  #<-- end `em()`
################################################################################















################################################################################
#' @rdname reml
#' @export
ai <- function(nuvin, skel, thetaG,
		   modMats, W, sLc, sln, r,
		   thetaR = NULL,  #<-- non-NULL if lambda==FALSE
		   sigma2e = NULL){  #<-- non-NULL if lambda==TRUE

  lambda <- is.null(thetaR)
  p <- length(nuvin)
  nuin <- vech2matlist(nuvin, skel)

  if(lambda){
    Rinv <- as(solve(matrix(sigma2e)), "symmetricMatrix")
  } else{
      Rinv <- as(solve(nuin[[thetaR]]), "symmetricMatrix")
    }
  B <- matrix(0, nrow = modMats$ny, ncol = p)

  # if no variance components (except residual), skip to the end
  if(length(thetaG) > 0){
    Ginv <- lapply(thetaG, FUN = function(x){as(solve(nuin[[x]]), "symmetricMatrix")})
    si <- modMats$nb+1
    for(g in thetaG){ #FIXME assumes thetaG is same length as nuvin
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # Meyer 1997: eqn. 20 [also f(theta) in Johnson & Thompson 1995 eqn 9c)]
      # Knight 2008: eqn 2.32-2.35 and 3.11-3.13 (p.34)
      #TODO for covariances see Johnson & Thompson 1995 eqn 11b
      Bg <- modMats$Zg[[g]] %*% sln[si:ei, , drop = FALSE] %*% Ginv[[g]] #FIXME is order correct? See difference in order between Johnson & Thompson (e.g., eqn. 11b) and Meyer 1997 eqn 20
      B[cbind(Bg@i+1, g)] <- Bg@x
      si <- ei+1
    }  #<-- end `for g`
  } else g <- p-1  #<-- end if no varcomps TODO: `g` assignment temporary until >1 residual


  #FIXME TODO Check what to do if more than 1 residual variance parameter
  if(g < p){
    if(lambda){
      B[, p] <- (modMats$y %*% Rinv)@x
    } else{
        B[, p] <- (r %*% Rinv)@x
      }
  }
  # Set up modified MME like the MMA of Meyer 1997 eqn. 23
  ## Substitute `B` instead of `y`
  if(lambda){
    BRHS <- Matrix(crossprod(W, B), sparse = TRUE) 
    tBRinvB <- crossprod(B)
  } else{
      # could pass KRinv and tWKRinv (if lambda=FALSE) from reml() into ai() 
      KRinv <- kronecker(Rinv, Diagonal(x = 1, n = modMats$Zr@Dim[[2L]]))
      tWKRinv <- crossprod(W, KRinv)
      BRHS <- Matrix(tWKRinv %*% B, sparse = TRUE)
      tBRinvB <- crossprod(B, KRinv) %*% B
    }
  ## tBPB
  ### Johnson & Thompson 1995 eqns 8,9b,9c
  #### (accomplishes same thing as Meyer 1997 eqns 22-23 for a Cholesky
  #### factorization as Boldman & Van Vleck eqn 6 applied to AI calculation)
  ## how to do this in c++ with Csparse
#   tS <- matrix(NA, nrow = ncol(BRHS), ncol = nrow(BRHS))
#    for(c in 1:p){
#      slv1P <- solve(sLc, BRHS[, c], system = "P")
#      slv1L <- solve(sLc, slv1P, system = "L")
#      slv1Lt <- solve(sLc, slv1L, system = "Lt")
#      tS[c, ] <- solve(sLc, slv1Lt, system = "Pt")@x
#    }
#    tBPB <- tBRinvB - (tS %*% BRHS)
  tBPB <- tBRinvB - crossprod(solve(sLc, BRHS, system = "A"), BRHS)
  AI <- 0.5 * tBPB 
  if(lambda) AI <- AI / sigma2e

 as(AI, "matrix")
}  #<-- end `ai()`
################################################################################    









################################################################################
#' @rdname reml
#' @param nb The number of columns in X.
#' @param listGeninv A list of generalized inverse matrices.
#' @param pinv An integer vector of the matrix permutation.
#' @export
tugug_trace <- function(thetaG, nb, rfxlvls, listGeninv, Cinv, sln, pinv = NULL){

  if(!is.null(pinv)){
    # do Cinv_ii before reorder Cinv
    Cinv_ii <- Cinv@x[attr(Cinv, "Zdiagp")][invPerm(pinv)] 
      P <- as(pinv, "pMatrix")
    Cinv <- crossprod(P, Cinv) %*% P
  } else Cinv_ii <- NULL
       
  # `trCinvGeninv_gg` = trace[Cinv_gg %*% geninv_gg]
  # `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
  ## `g` is the gth component of the G-structure to model
  ## `geninv` is the generalized inverse
  ### (not the inverse of the G-matrix/(co)variance components)
  trCinvGeninv_gg <- tugug <- as.list(rep(0, length(thetaG)))
  si <- nb + 1
  for(g in thetaG){
    qi <- rfxlvls[g]
    ei <- si - 1 + qi
#TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    ## Johnson & Thompson 1995 equations 9 & 10
    if(inherits(listGeninv[[g]], "ddiMatrix")){
    ### No generalized inverse associated with the random effects
    ### Johnson & Thompson 1995 eqn 10a
      #### 3rd term in the equation
      tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE])
      #### 2nd term in the equation
      trCinvGeninv_gg[[g]] <- ifelse(is.null(Cinv_ii),
          tr(Cinv[si:ei, si:ei]), sum(Cinv_ii[si:ei]))
     
    } else{
      ### Generalized inverse associated with the random effects
      ### Johnson & Thompson 1995 eqn 9a
        #### 3rd term in the equation
        tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE],
                         listGeninv[[g]]) %*% sln[si:ei, , drop = FALSE]
        #### 2nd term in the equation
        # the trace of a product = the sum of the element-by-element product
        ## trace(geninv[[g]] %*% Cinv[si:ei, si:ei]
        ###                       = sum((geninv[[g]] * Cinv[si:ei, si:ei])@x)
        trCinvGeninv_gg[[g]] <- sum((listGeninv[[g]] * Cinv[si:ei, si:ei])@x)
      }  #<-- end if else whether a diagonal g-inverse associated with rfx

    si <- ei+1

  }  #<-- end `for g in thetaG`

  return(list(tugug = tugug,
              trace = trCinvGeninv_gg))
}  #<-- end tugug_trace              






################################################################################
#' @rdname reml
#' @export
gradFun <- function(nuvin, thetaG, rfxlvls, sln,
	tugug, trace,
	sigma2e = NULL,   #<-- non-NULL if lambda==TRUE
	r = NULL, nminfrfx = NULL){  #<-- non-NULL if lambda==FALSE

  lambda <- is.null(r)  #<-- lambda scale model
  p <- length(nuvin)
  dLdnu <- matrix(0, nrow = p, ncol = 1, dimnames = list(names(nuvin), NULL))
  # tee = e'e
  if(!lambda) tee <- crossprod(r)

  # Skip most of below if there are no variance components other than residual
  if(length(thetaG) > 0){

    # First derivatives (gradient/score)
    if(lambda){
      for(g in thetaG){
        # Johnson and Thompson 1995 Appendix 2 eqn B3 and eqn 9a and 10a
        dLdnu[g] <- (rfxlvls[g] / nuvin[g]) - (1 / nuvin[g]^2) * (trace[[g]] + tugug[[g]] / sigma2e)
      }  #<-- end for g

    } else{  #<-- else when NOT lambda scale
        ## Johnson and Thompson 1995 eqn 9b
#FIXME change `[p]` below to be number of residual (co)variances
        dLdnu[p] <- (nminfrfx / tail(nuvin, 1)) - (tee / tail(nuvin, 1)^2)
        for(g in thetaG){
          dLdnu[p] <- dLdnu[p] + (1 / tail(nuvin, 1)) * (trace[[g]] /nuvin[g]) 
          # Johnson and Thompson 1995 eqn 9a and 10a
          dLdnu[g] <- (rfxlvls[g] / nuvin[g]) - (1 / nuvin[g]^2) * (trace[[g]] + tugug[[g]])
        }  #<-- end for g
      }
  } else{  #<-- end if varcomps besides residuals
      # First derivatives (gradient/score) if only residual variance in model
      ## if lambda==TRUE then don't do this at all and `dLdnu` is empty
###FIXME when >1 residual (co)variance need to change whether skip with lambda=TRUE
      if(!lambda) dLdnu[p] <- (nminfrfx / tail(nuvin, 1)) - (tee / tail(nuvin, 1)^2)
    }  #<-- end if/else no varcomps besides residual

 # Johnson and Thompson 1995 don't use -0.5, because likelihood is -2 log likelihood
 ## see `-2` on left-hand side of Johnson & Thompson eqn 3
 -0.5 * dLdnu
}  #<-- end `gradFun()`
##########################
# Changing ginverse elements to `dsCMatrix` doesn't speedup traces, since
## they end up more or less as dense matrices but in dgCMatrix from the product
##########################



################################################################################
#' @rdname reml
#' @param grObj An list of class \code{grMod}.
#' @param lL A numeric value for REML log-likelihood value.
#' @param fd A character indicating whether forward, combined, or backward finite
#'   differences (\dQuote{fdiff}, \dQuote{cdiff}, or \dQuote{bdiff}, respectively)
#'   are to be calculated.
#' @export
gradFun_fd <- function(nuvin, grObj, lL, fd = c("fdiff", "cdiff", "bdiff")){

  fd <- match.arg(fd)
  h <- grObj$h
  denomSC <- ifelse(fd == "cdiff", 2, 1)
    
  nuv_tmp <- nuvin
    skel <- attr(nuvin, "skel")
  lambda <- grObj$lambda
    if(lambda){
      thetaR <- NULL
      tWW <- grObj$tWW
      RHS <- grObj$RHS
    } else{
        thetaR <- grObj$thetaR
        tWW <- RHS <- NULL
      }
  thetaG <- grObj$thetaG
  p <- grObj$p
  con <- grObj$con
  fxL <- fxU <- matrix(lL, nrow = grObj$p, ncol = 1,
    dimnames = list(names(nuvin), NULL))

  # Skip most of below if there are no variance components other than residual
  if(length(thetaG) > 0){
    # `g` is the gth component of the G-structure to model
    for(g in thetaG){
      if(con[g] != "F"){
        if(fd != "bdiff"){  # if either FORWARD or CENTRAL difference method
          nuv_tmp[g] <- nuvin[g] + h
          lL_fd <- reml(vech2matlist(nuv_tmp, skel), skel, thetaG,
	    grObj$sLc, grObj$modMats, grObj$W, grObj$Bpinv,
            grObj$nminffx, grObj$nminfrfx, grObj$rfxlvls,
            grObj$rfxIncContrib2loglik,
	    thetaR, tWW, RHS)$loglik
	  fxL[g, 1] <- lL_fd   # either forward or central diff.
        }     
        if(fd != "fdiff"){  # if either BACKWARD or CENTRAL difference method
          nuv_tmp[g] <- nuvin[g] - h
          lL_fd <- reml(vech2matlist(nuv_tmp, skel), skel, thetaG,
	    grObj$sLc, grObj$modMats, grObj$W, grObj$Bpinv,
            grObj$nminffx, grObj$nminfrfx, grObj$rfxlvls,
            grObj$rfxIncContrib2loglik,
	    thetaR, tWW, RHS)$loglik
	  fxU[g, 1] <- lL_fd   # either backward or central diff.
        }
        
        nuv_tmp[g] <- nuvin[g]  #<-- reset
      }  #<-- end if g is not fixed
    }  #<-- end `for g in thetaG`
  }  #<-- end if there are any G-structure elements


  # Residual (co)variances when not on Lambda scale
#FIXME change `[p]` below to be number of residual (co)variances
  if(!lambda){
    if(con[p] != "F"){
      if(fd != "bdiff"){  # if either FORWARD or CENTRAL difference method
        nuv_tmp[p] <- nuvin[p] + h
        lL_fd <- reml(vech2matlist(nuv_tmp, skel), skel, thetaG,
	  grObj$sLc, grObj$modMats, grObj$W, grObj$Bpinv,
          grObj$nminffx, grObj$nminfrfx, grObj$rfxlvls,
          grObj$rfxIncContrib2loglik,
	  thetaR, tWW, RHS)$loglik
	# residual is backward hence `fxL` (instead of `fxU`)  
        fxL[p, 1] <- lL_fd   # either forward or central diff.
      }     
      if(fd != "fdiff"){  # if either BACKWARD or CENTRAL difference method
        nuv_tmp[p] <- nuvin[p] - h
        lL_fd <- reml(vech2matlist(nuv_tmp, skel), skel, thetaG,
	  grObj$sLc, grObj$modMats, grObj$W, grObj$Bpinv,
          grObj$nminffx, grObj$nminfrfx, grObj$rfxlvls,
          grObj$rfxIncContrib2loglik,
	  thetaR, tWW, RHS)$loglik
	# residual is backward hence `fxU` (instead of `fxL`) 
        fxU[p, 1] <- lL_fd   # either backward or central diff.
      }
      
      nuv_tmp[p] <- nuvin[p]  #<-- reset
    }  #<-- end if g is not fixed
  }  #<-- end if lambda=FALSE    

  # First derivatives (gradient/score)
  ## forward = [f(x+h) - f(x)] / h
  ## backward = [f(x) - f(x-h)] / h
  ## central = [f(x+h) - f(x-h)] / 2h
  dLdnu <- (fxU - fxL) / (denomSC * h)

 -1 * dLdnu  #<-- since optimizing the negative log-likelihood
}  #<-- end `gradFun_fd()`

