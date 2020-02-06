###############################
# reml_lambda()
###############################
#' @rdname gremlinRmod
#' @export
reml_lambda <- function(nu, skel, thetaG, sLc,
	modMats, W, tWW, Bpinv, RHS,
	nminffx, nminfrfx, rfxlvls, rfxIncContrib2loglik){

  ##1c Now make coefficient matrix of MME
  ##`sapply()` to invert G_i and multiply with ginverse element (e.g., I or geninv)
  if(modMats$nG > 0){
    C <- as(tWW + bdiag(c(Bpinv,
      sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGeninv[[u]], solve(nu[[u]]))}))), "symmetricMatrix")
    } else C <- as(tWW + Bpinv, "symmetricMatrix")
  ## Update symbolic factorization with variance parameters from this iteration
  sLc <- update(sLc, C)

  # 5 record log-like, check convergence, & determine next varcomps to evaluate  
  ##5a determine log(|C|) and y'Py
  # Boldman and Van Vleck 1991 eqn. 6 (tyPy) and 7 (log|C|)    
  # Meyer 1997 eqn. 13 (tyPy)
  tyPy <- crossprod(modMats$y) - crossprod(solve(sLc, RHS, system = "A"), RHS)
  logDetC <- 2 * sum(log(sLc@x[sLc@p+1][1:sLc@Dim[[1L]]]))
  # alternatively see `determinant` method for CHMfactor

  # Residual variance
  sigma2e <- tyPy / nminffx

  # Construct the log-likelihood (Meyer 1997, eqn. 8)
  ## (first put together as `-2log-likelihood`)
  ## 'log(|R|)'
  ### assumes X of full rank
  loglik <- nminfrfx * log(sigma2e)

  # 'log(|G|)'
  #FIXME: Only works for independent random effects right now!
  logDetGfun <- function(x){rfxlvls[x] * log(as.vector(nu[[x]]*sigma2e))}
  loglik <- loglik + if(modMats$nG == 0){ 0
    } else{
        sum(sapply(seq(modMats$nG), FUN = logDetGfun)) + rfxIncContrib2loglik
      }

  ## log(|C|) + tyPy
  ### nminffx is tyPy/sigma2e, simplified because sigma2e = tyPy / nminffx
  ### and mulitply by -0.5 to calculate `loglik` from `-2loglik`
  loglik <- -0.5 * (loglik + logDetC + nminffx)




  # solve MME for BLUEs/BLUPs
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

  # calculate residuals
  r <- modMats$y - W %*% sln


 return(structure(list(loglik = loglik@x,
		sigma2e = sigma2e@x, tyPy = tyPy@x, logDetC = logDetC,
		sln = sln, r = r, sLc = sLc),
	class = "gremlin"))
}  #<-- end `reml_lambda()` 


################################################################################
















###############################
# reml() 
###############################
#' @rdname gremlinRmod
#' @export
reml <- function(nu, skel, thetaG, thetaR, sLc,
	modMats, W, Bpinv, nminffx, nminfrfx, rfxlvls, rfxIncContrib2loglik){
  Rinv <- as(solve(nu[[thetaR]]), "symmetricMatrix")
  Ginv <- lapply(thetaG, FUN = function(x){as(solve(nu[[x]]), "symmetricMatrix")}) # Meyer 1991, p.77 (<-- also see MMA on p70/eqn4)

  
  ##1c Now make coefficient matrix of MME
  ### Rinv Kronecker with Diagonal/I
  KRinv <- kronecker(Rinv, Diagonal(x = 1, n = modMats$Zr@Dim[[2L]]))
  tyRinvy <- as(crossprod(modMats$y, KRinv) %*% modMats$y, "sparseMatrix")
    #TODO: what to do if y is multivariate (still just 1 column, so D just 1 number?
    ### Components of Meyer 1989 eqn 2
    tWKRinv <- crossprod(W, KRinv)
    tWKRinvW <- tWKRinv %*% W  
  ##`sapply()` to multiply inverse G_i with ginverse element (e.g., I or geninv)
  if(modMats$nG > 0){
    C <- as(tWKRinvW + bdiag(c(Bpinv,
      sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGeninv[[u]], Ginv[[u]])}))), "symmetricMatrix")
  } else C <- as(tWKRinvW + Bpinv, "symmetricMatrix")
  ## Update symbolic factorization with variance parameters from this iteration
  sLc <- update(sLc, C)
  RHS <- Matrix(tWKRinv %*% modMats$y, sparse = TRUE)

  # 5 record log-like, check convergence, & determine next varcomps to evaluate  
  ##5a determine log(|C|) and y'Py
  # Boldman and Van Vleck 1991 eqn. 6 (tyPy) and 7 (log|C|)    
  # Meyer 1997 eqn. 13 (tyPy)
  tyPy <- tyRinvy - crossprod(solve(sLc, RHS, system = "A"), RHS)
  logDetC <- 2 * sum(log(sLc@x[sLc@p+1][1:sLc@Dim[[1L]]]))
  # alternatively see `determinant` method for CHMfactor

  # Construct the log-likelihood (Meyer 1997, eqn. 8)
  ## (first put together as `-2log-likelihood`)
  ## 'log(|R|)' Meyer 1997 eqn. 10
  ### assumes X of full rank
  loglik <- modMats$ny * log(nu[[thetaR]])

  ## 'log(|G|)' Meyer 1997 eqn. 9
  #FIXME: Only works for independent random effects right now!
  logDetGfun <- function(x){rfxlvls[x] * log(as.vector(nu[[x]]))}
  loglik <- loglik + if(modMats$nG == 0){ 0
    } else{
        sum(sapply(seq(modMats$nG), FUN = logDetGfun)) + rfxIncContrib2loglik
      }

  ## log(|C|) + tyPy
  ### and mulitply by -0.5 to calculate `loglik` from `-2loglik`
  loglik <- -0.5 * (loglik + logDetC + tyPy)




  # solve MME for BLUEs/BLUPs
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

  # calculate residuals
  r <- modMats$y - W %*% sln

 return(structure(list(loglik = loglik@x,
		sigma2e = 1, tyPy = tyPy@x, logDetC = logDetC,
		sln = sln, r = r, sLc = sLc),
	class = "gremlin"))
}  #<-- end `reml()` 


################################################################################



















################################################################################
## Meyer and Smith 1996 for algorithm using derivatives of loglik
### eqn 15-18 (+ eqn 33-42ish) for derivatives of tyPy and logDetC
### Smith 1995 for very technical details
# EM refs: Hofer 1998 eqn 10-12
## XXX note Hofer eqn 12 missing sigma2e in last term of non-residual formula
### see instead Mrode 2005 (p. 241-245)
#' @rdname gremlinRmod
#' @export
em <- function(nuvin, thetaG, thetaR,
	modMats, nminffx, sLc, ndgeninv, sln, r){
  ## go "backwards" so can fill in Lc with lower triangle of Cinv
  ei <- modMats$nb + sum(sapply(modMats$Zg, FUN = ncol))
  Ig <- Diagonal(n = sLc@Dim[1L], x = 1)
  Cinv_ii <- matrix(0, nrow = nrow(sln), ncol = 1)
  for(g in rev(thetaG)){
    qi <- ncol(modMats$Zg[[g]])
    si <- ei - qi + 1
    # Note: trace of a product == the sum of the element-by-element product
    ## Consequently, don't have to make `Cinv`, just diagonals
#XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
##XXX eqn. 2.44 is the score/gradient! for a varcomp
    trace <- 0
    if(ndgeninv[g]){
      o <- (crossprod(sln[si:ei, , drop = FALSE], modMats$listGeninv[[g]]) %*% sln[si:ei, , drop = FALSE])@x
      for(k in si:ei){
        Cinv_siei_k <- solve(sLc, b = Ig[, k], system = "A")[si:ei, , drop = TRUE]
        Cinv_ii[k] <- Cinv_siei_k[k-si+1]       
        trace <- trace + sum(modMats$listGeninv[[g]][(k-si+1), , drop = TRUE] * Cinv_siei_k)
      }  #<-- end for k
    } else{
        ## first term
        o <- crossprod(sln[si:ei, , drop = FALSE])
        for(k in si:ei){
          Cinv_ii[k] <- solve(sLc, b = Ig[, k], system = "A")[k,]
          trace <- trace + Cinv_ii[k]
        }  #<-- end for k
      }  #<-- end if/else ndgeninv
    #TODO check `*tail(nuv,1)` correctly handles models with covariance matrices
browser()
    nuvin[g] <- as(as(matrix((o + trace) / qi),  #XXX was `trace*tail(nuvin,1)` 
      "symmetricMatrix"),
      "dsCMatrix")
    ei <- si-1
  }  #<-- end `for g`

  nuvin[thetaR] <- crossprod(modMats$y, r) / nminffx

 return(structure(list(nuv = nuvin, Cinv_ii = Cinv_ii),
	class = "gremlin"))
}  #<-- end `em()`
################################################################################















################################################################################
#' @rdname gremlinRmod
#' @export
#Meyer 1996:
## Likelihood eqn 3 <--> log determinants
#XXX `ai_lambda()` used for when Rinv HAS BEEN FACTORED OUT of MME
ai_lambda <- function(nuvin, skel, thetaG, sigma2e,
		   modMats, W, sLc, sln, r){
  p <- length(nuvin)
  nuin <- vech2matlist(nuvin, skel)
  Rinv <- as(solve(matrix(sigma2e)), "symmetricMatrix")
  Ginv <- lapply(thetaG, FUN = function(x){as(solve(nuin[[x]]), "symmetricMatrix")})
  si <- modMats$nb+1
  B <- Matrix(0, nrow = modMats$ny, ncol = p, sparse = TRUE)
  for(g in thetaG){ #FIXME assumes thetaG is same length as nuvin
    qi <- ncol(modMats$Zg[[g]])
    ei <- si - 1 + qi
    # Meyer 1997: eqn. 20 [also f(theta) in Johnson & Thompson 1995 eqn 9c)]
    # Knight 2008: eqn 2.32-2.35 and 3.11-3.13 (p.34)
    #TODO for covariances see Johnson & Thompson 1995 eqn 11b
    B[, g] <- modMats$Zg[[g]] %*% sln[si:ei, , drop = FALSE] %*% Ginv[[g]] #FIXME is order correct? See difference in order between Johnson & Thompson (e.g., eqn. 11b) and Meyer 1997 eqn 20
    si <- ei+1
  }  #<-- end `for g`

  #FIXME TODO Check what to do if more than 1 residual variance parameter
  if(g < p){
    B[, p] <- modMats$y %*% Rinv
  }
  # Set up modified MME like the MMA of Meyer 1997 eqn. 23
  ## Substitute `B` instead of `y`
  BRHS <- Matrix(crossprod(W, B), sparse = TRUE)
  ## tBPB
  ### Meyer 1997 eqns 22-23 (extends Johnson & Thompson 1995 eqns 8,9b,9c)
  tBPB <- crossprod(B) - crossprod(solve(sLc, BRHS, system = "A"), BRHS)
  AI <- 0.5 * tBPB / sigma2e

 return(structure(list(AI = as(AI, "matrix")),
	class = "gremlin"))
}  #<-- end `ai_lambda()`
################################################################################    













################################################################################
#' @rdname gremlinRmod
#' @export
#Meyer 1996:
## Likelihood eqn 3 <--> log determinants
#XXX `ai()` used for when Rinv is NOT factored out of MME
ai <- function(nuvin, skel, thetaG, thetaR,
		   modMats, W, sLc, sln, r){
  p <- length(nuvin)
  nuin <- vech2matlist(nuvin, skel)
  Rinv <- as(solve(nuin[[thetaR]]), "symmetricMatrix")
  Ginv <- lapply(thetaG, FUN = function(x){as(solve(nuin[[x]]), "symmetricMatrix")})
  si <- modMats$nb+1
  B <- Matrix(0, nrow = modMats$ny, ncol = p, sparse = TRUE)
  for(g in thetaG){ #FIXME assumes thetaG is same length as nuvin
    qi <- ncol(modMats$Zg[[g]])
    ei <- si - 1 + qi
    # Meyer 1997: eqn. 20 [also f(theta) in Johnson & Thompson 1995 eqn 9c)]
    # Knight 2008: eqn 2.32-2.35 and 3.11-3.13 (p.34)
    #TODO for covariances see Johnson & Thompson 1995 eqn 11b
    B[, g] <- modMats$Zg[[g]] %*% sln[si:ei, , drop = FALSE] %*% Ginv[[g]] #FIXME is order correct? See difference in order between Johnson & Thompson (e.g., eqn. 11b) and Meyer 1997 eqn 20
    si <- ei+1
  }  #<-- end `for g`

  #FIXME TODO Check what to do if more than 1 residual variance parameter
  if(g < p){
    B[, p] <- r %*% Rinv
  }
  # Set up modified MME like the MMA of Meyer 1997 eqn. 23
  ## Substitute `B` instead of `y`
  KRinv <- kronecker(Rinv, Diagonal(x = 1, n = modMats$Zr@Dim[[2L]]))
  tWKRinv <- crossprod(W, KRinv)
  tBRinvB <- crossprod(B, KRinv) %*% B
  BRHS <- Matrix(tWKRinv %*% B, sparse = TRUE)
  ## tBPB
  ### Meyer 1997 eqns 22-23 (extends Johnson & Thompson 1995 eqns 8,9b,9c)
  tBPB <- tBRinvB - crossprod(solve(sLc, BRHS, system = "A"), BRHS)
  AI <- 0.5 * tBPB 

 return(structure(list(AI = as(AI, "matrix")),
	class = "gremlin"))
}  #<-- end `ai()`
################################################################################    








################################################################################
#' @rdname gremlinRmod
#' @export
gradFun_lambda <- function(nuvin, thetaG, modMats, Cinv, sln, sigma2e){
  p <- length(nuvin)
  dLdnu <- matrix(0.0, nrow = p, ncol = 1, dimnames = list(names(nuvin), NULL))

  # `trCinvGeninv_gg` = trace[Cinv_gg %*% geninv_gg]
  # `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
  ## `g` is the gth component of the G-structure to model
  ## `geninv` is the generalized inverse
  ### (not the inverse of the G-matrix/(co)variance components)
  trCinvGeninv_gg <- tugug <- as.list(rep(0, length(thetaG)))
  si <- modMats$nb+1
  for(g in thetaG){
    qi <- ncol(modMats$Zg[[g]])
    ei <- si - 1 + qi
#TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    ## Johnson & Thompson 1995 equations 9 & 10
    if(class(modMats$listGeninv[[g]]) == "ddiMatrix"){
    ### No generalized inverse associated with the random effects
    ### Johnson & Thompson 1995 eqn 10a
      #### 3rd term in the equation
      tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE])
      #### 2nd term in the equation
      trCinvGeninv_gg[[g]] <- tr(Cinv[si:ei, si:ei])
       
    } else{
      ### Generalized inverse associated with the random effects
      ### Johnson & Thompson 1995 eqn 9a
        #### 3rd term in the equation
        tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE], modMats$listGeninv[[g]]) %*% sln[si:ei, , drop = FALSE]
        #### 2nd term in the equation
#TODO: DOES the trace of a product = the sum of the element-by-element product?
        trCinvGeninv_gg[[g]] <- tr(modMats$listGeninv[[g]] %*% Cinv[si:ei, si:ei])
      }  #<-- end if else whether a diagonal g-inverse associated with rfx

    si <- ei+1
  }  #<-- end `for g in thetaG`

  # First derivatives (gradient/score)
  for(g in thetaG){
    # Johnson and Thompson 1995 Appendix 2 eqn B3 and eqn 9a and 10a
    dLdnu[g] <- (ncol(modMats$Zg[[g]]) / nuvin[g]) - (1 / nuvin[g]^2) * (trCinvGeninv_gg[[g]] + tugug[[g]] / sigma2e)
  }
 # Johnson and Thompson 1995 don't use -0.5, because likelihood is -2 log likelihood
 ## see `-2` on left-hand side of Johnson & Thompson eqn 3
 -0.5 * dLdnu
}  #<-- end `gradFun_lambda()`
##########################
# Changing ginverse elements to `dsCMatrix` doesn't speedup traces, since
## they end up more or less as dense matrices but in dgCMatrix from the product
##########################








################################################################################
#' @rdname gremlinRmod
#' @export
gradFun <- function(nuvin, thetaG, modMats, Cinv, nminfrfx, sln, r){
  p <- length(nuvin)
  dLdnu <- matrix(NA, nrow = p, ncol = 1, dimnames = list(names(nuvin), NULL))
  # tee = e'e
  tee <- crossprod(r)

  # `trCinvGeninv_gg` = trace[Cinv_gg %*% geninv_gg]
  # `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
  ## `g` is the gth component of the G-structure to model
  ## `geninv` is the generalized inverse
  ### (not the inverse of the G-matrix/(co)variance components)
  trCinvGeninv_gg <- tugug <- as.list(rep(0, length(thetaG)))
  si <- modMats$nb+1
  for(g in thetaG){
    qi <- ncol(modMats$Zg[[g]])
    ei <- si - 1 + qi
#TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    ## Johnson & Thompson 1995 equations 9 & 10
    if(class(modMats$listGeninv[[g]]) == "ddiMatrix"){
    ### No generalized inverse associated with the random effects
    ### Johnson & Thompson 1995 eqn 10a
      #### 3rd term in the equation
      tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE])
      #### 2nd term in the equation
      trCinvGeninv_gg[[g]] <- tr(Cinv[si:ei, si:ei])
      
    } else{
      ### Generalized inverse associated with the random effects
      ### Johnson & Thompson 1995 eqn 9a
        #### 3rd term in the equation
        tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE], modMats$listGeninv[[g]]) %*% sln[si:ei, , drop = FALSE]
        #### 2nd term in the equation
#TODO: DOES the trace of a product = the sum of the element-by-element product?
        trCinvGeninv_gg[[g]] <- tr(modMats$listGeninv[[g]] %*% Cinv[si:ei, si:ei])
      }  #<-- end if else whether a diagonal g-inverse associated with rfx

    si <- ei+1
  }  #<-- end `for g in thetaG`

  # First derivatives (gradient/score)
#FIXME change `[p]` below to be number of residual (co)variances
  ## Johnson and Thompson 1995 eqn 9b
  dLdnu[p] <- (nminfrfx / tail(nuvin, 1)) - (tee / tail(nuvin, 1)^2)
  for(g in thetaG){
    dLdnu[p] <- dLdnu[p] + (1 / tail(nuvin, 1)) * (trCinvGeninv_gg[[g]] /nuvin[g]) 
    # Johnson and Thompson 1995 eqn 9a and 10a
    dLdnu[g] <- (ncol(modMats$Zg[[g]]) / nuvin[g]) - (1 / nuvin[g]^2) * (trCinvGeninv_gg[[g]] + tugug[[g]])
  }
 # Johnson and Thompson 1995 don't use -0.5, because likelihood is -2 log likelihood
 ## see `-2` on left-hand side of Johnson & Thompson eqn 3
 -0.5 * dLdnu
}  #<-- end `gradFun()`
##########################
# Changing ginverse elements to `dsCMatrix` doesn't speedup traces, since
## they end up more or less as dense matrices but in dgCMatrix from the product
##########################

