###############################
# reml()
###############################
#' @rdname gremlinRmod
#' @export
reml <- function(thetav, skel, thetaG, thetaR,
	modMats, W, tWW, Bpinv, RHS, D,
	nminffx, nminfrfx, rfxlvls, rfxIncContrib2loglik, rm1, sLc){
      theta <- vech2matlist(thetav, skel)
      nu <- theta
#XXX Don't need for EM algorithm
#    transform <- FALSE
#    if(transform){
#      nu <- list(G = lapply(theta$G, FUN = function(x){log(chol(x))}),
#		R = log(chol(theta$R)))
#    } else nu <- theta
      cRinv <- solve(chol(nu[[thetaR]]))
      Ginv <- lapply(thetaG, FUN = function(x){as(crossprod(cRinv, nu[[x]]) %*% cRinv, "symmetricMatrix")}) # Meyer 1991, p.77 (<-- ?p70/eqn4?)
  

      ##1c Now make coefficient matrix of MME
      ##`sapply()` to invert G_i and multiply with ginverse element (e.g., I or geninv)
      if(modMats$nG > 0){
        C <- as(tWW + bdiag(c(Bpinv,
	  sapply(1:modMats$nG, FUN = function(u){kronecker(modMats$listGeninv[[u]], solve(Ginv[[u]]))}))), "symmetricMatrix")
      } else C <- as(tWW + Bpinv, "symmetricMatrix")

      M <- as(cbind(rbind(C, t(RHS)),
	    rbind(RHS, D)), "symmetricMatrix")





#FIXME XXX TURN Permutation back ON perm==TRUE
      sLm <- Cholesky(M, perm = FALSE, LDL = FALSE, super = FALSE)

     # 5 record log-like, check convergence, & determine next varcomps to evaluate  
      ##5a determine log(|C|) and y'Py
      ### Meyer & Smith 1996, eqns 12-14 (and 9)
      #### Also see Meyer & Kirkpatrick 2005 GSE. eqn. 18: if cholesky of MMA = LL'
      # Meyer & Smith 1996, eqn. 14
      tyPy <- tail(sLm@x, 1)^2
      logDetC <- 2 * sum(log(sLm@x[sLm@p+1][1:(sLm@Dim[[1L]]-1)]))
      # alternatively from sLc: `... sum(log(sLc@x[sLc@p+1][1:sLc@Dim[[1]]]))`
      ## XXX does sLm takes longer than sLc for more complex models

      # Residual variance
      sigma2e <- tyPy / nminffx

      # nminffx is tyPy/sigma2e, simplified because sigma2e = tyPy / nminffx
      loglik <- nminffx + logDetC

      # 'log(|R|)'
      #TODO: assumes X of full rank
      loglik <- loglik + nminfrfx * log(sigma2e)
      #ALTERNATIVELY: If Rinv NOT factored out of MMA `loglik <- loglik + ny * log(Rstart)`

      # 'log(|G|)'
      #FIXME: Only works for independent random effects right now!
      loglik <- -0.5 * (loglik + if(modMats$nG == 0) 0 else sum(sapply(seq(modMats$nG), FUN = function(x){rfxlvls[x] * log(as.vector(Ginv[[x]]*sigma2e))})) + rfxIncContrib2loglik)
      # Below uses original starting value for residual variances - for agreement with WOMBAT
      #loglik <- -0.5 * (loglik + sum(sapply(seq(nG), FUN = function(x){rfxlvls[x] * log(as.vector(Gstart[[x]]))})) + rfxIncContrib2loglik)






      # solve MME for BLUEs/BLUPs
#XXX Do I need to solve for BLUEs (see Knight 2008 for just getting BLUPs eqn 2.13-15)
      ## see Mrode 2005 chapter
      # permuted RHS of MME is `RHSperm`
      #sln[] <<- Pc %*% solve(a = sLc, b = t(RHSperm), system = "A")
      #otherwise
      #sln[] <<- solve(a = C, b = RHS, system = "A")
      # alternatively, get factorization of C (sLc) from sLm:
      ## Assume same non-zero pattern
      sLc@x <- sLm@x[rm1]

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

     return(structure(list(loglik = loglik,
		sigma2e = sigma2e, tyPy = tyPy, logDetC = logDetC,
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
  em <- function(thetain, thetaG, thetaR,
	modMats, nminffx, sLc, ndgeninv, sln, r){
    ## go "backwards" so can fill in Lc with lower triangle of Cinv
    ei <- modMats$nb + sum(sapply(modMats$Zg, FUN = ncol))
    Ig <- Diagonal(n = sLc@Dim[1L], x = 1)
    Cinv_ii <- matrix(0, nrow = nrow(sln), ncol = 1)
    for(g in rev(thetaG)){
      qi <- ncol(modMats$Zg[[g]])
      si <- ei - qi + 1
      # note the trace of a product is equal to the sum of the element-by-element product
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
      thetav <- sapply(thetain, FUN = slot, name = "x")
      #TODO check that this correctly handles models with covariance matrices
      thetain[[g]] <- as(as(matrix((o + trace*tail(thetav, 1)) / qi), "symmetricMatrix"), "dsCMatrix")
      ei <- si-1
    }
    #FIXME make sure `nminffx` == `ncol(X)` even when reduced rank
    thetain[[thetaR]] <- crossprod(modMats$y, r) / nminffx
   return(structure(list(theta = thetain, Cinv_ii = Cinv_ii),
	class = "gremlin"))
  }
################################################################################















################################################################################
#' @rdname gremlinRmod
#' @export
  ai <- function(thetavin, skel, thetaG, thetaR, modMats, W, Cinv, nminfrfx,
	sln, r, ezero){
    p <- length(thetavin)
    thetain <- vech2matlist(thetavin, skel)
    # setup Vinv directly [See Johnson & Thompson 1995, Appendix 1 (eqn A1)]
##FIXME will it work for >1 residual (co)variance????  XXX

    #TODO below removes X from W, make sure `nb` is correct when reduced rank X
    tmptmpGinv <- with(modMats, bdiag(sapply(seq(nG), FUN = function(u){kronecker(listGeninv[[u]], solve(thetain[[u]]))})))
    #FIXME why the difference between using thetain versus Ginv
    #tmptmpGinv <- with(modMats, bdiag(sapply(nG, FUN = function(u){kronecker(listGeninv[[u]], solve(Ginv[[u]]))})))
    Rinv <- kronecker(Diagonal(n = modMats$ny, x = 1), thetaR) #FIXME change thetaR to nu or transformed scale? 
    tZRinvZ <- with(modMats, crossprod(W[, -c(1:nb)], Rinv) %*%  W[, -c(1:nb)])
    #TODO can I use tWW below?
    PorVinv <- with(modMats, Rinv - tcrossprod(Rinv %*% W[, -c(1:nb)] %*% solve(tZRinvZ + tmptmpGinv), W[, -c(1:nb)]) %*% Rinv)  #<-- FIXME move outside?
#FIXME why P and P2 different?
    PorVinv <- with(modMats, PorVinv - PorVinv %*% X %*% tcrossprod(solve(crossprod(X, PorVinv) %*% X), X) %*% PorVinv)
#    P2 <- Rinv - tcrossprod(Rinv %*% W %*% Cinv, W) %*% Rinv #<-- See Gilmour et al. 1995 end of p.1441
#    P3 <- Diagonal(n = nrow(Rinv), x = 1) - W %*% solve(C) %*% t(W) #<-- AIreml_heritabilityPkg

    # tee = e'e
    tee <- crossprod(r)
    # trCinvGeninv_gg = trace[Cinv_gg %*% geninv_gg] | tugug = u_gg' %*% geninv_gg %*% u_gg
    ## g is the gth component of the G-structure to model
    ## geninv is the generalized inverse (not the inverse of the G-matrix/varcomps)
#TODO make variable `length(thetaG)`
    trCinvGeninv_gg <- tugug <- as.list(rep(0, length(thetaG)))
    si <- modMats$nb+1
    AI <- matrix(NA, nrow = p, ncol = p)
    dLdtheta <- matrix(NA, nrow = p, ncol = 1, dimnames = list(names(thetavin), NULL))
    for(g in thetaG){ #FIXME assumes thetaG is same length as thetavin
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      isln <- sln[si:ei, , drop = FALSE]
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      tugug[[g]] <- crossprod(isln, modMats$listGeninv[[g]]) %*% isln
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      trCinvGeninv_gg[[g]] <- tr(modMats$listGeninv[[g]] %*% Cinv[si:ei, si:ei])
#      A <- solve(Ainv)
#      AI[g, g] <- 0.5 * (t(y) %*% P %*% A %*% P %*% A %*% P %*% y) / thetavin[g]
#FIXME Check for multivariate when theta is a matrix, but below g is assumed to be a single (co)variance
      AI[g, g] <- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, si:ei] %*% isln) / (thetavin[g]^2))@x
#(crossprod(isln, crossprod(modMats$Zg[[g]], P)) %*% modMats$Zg[[g]] %*% isln) / (thetavin[g]^2)
      if((g+1) < p){
        for(k in (g+1):(p-1)){  #<-- fill upper triangle
          sk <- sum(sapply(seq(k-1), FUN = function(u){ncol(modMats$Zg[[u]])})) + modMats$nb + 1
          ek <- sk - 1 + ncol(modMats$Zg[[k]])
          AI[g, k] <- AI[k, g] <- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, sk:ek] %*% sln[sk:ek, , drop = FALSE]) / (thetavin[g]*thetavin[k]))@x    
        }  #<-- end 'for(k ...)`
      }  #<-- end `if()`
      AI[g, p] <- AI[p, g] <- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% r) / (thetavin[g]*thetavin[p]))@x
      si <- ei+1
    }   #<-- end `for(g ...)`
    AI[p, p] <- 0.5 * ((crossprod(r, PorVinv) %*% r) / thetavin[p]^2)@x

    # First derivatives (gradient)
#TODO check that `nminfrfx` is still `n-p-q` when more than 1 random effect 
## ALSO check what `q` of `n-p-q` is when >1 random effect
#FIXME change `[p]` below to be number of residual (co)variances
    dLdtheta[p] <- 0.5*((tee / thetavin[p]^2) - (nminfrfx) / thetavin[p]) 
    for(g in thetaG){
      dLdtheta[p] <- dLdtheta[p] - 0.5*(trCinvGeninv_gg[[g]]/thetavin[g]) 
      dLdtheta[g] <- 0.5*(tugug[[g]]/(thetavin[g]^2) - ncol(modMats$Zg[[g]])/thetavin[g] + trCinvGeninv_gg[[g]]*thetavin[p]/(thetavin[g]^2))
    }

    rcondAI <- rcond(AI)
    if(rcondAI < ezero){
      if(v > 1){
        cat("Reciprocal condition number of AI matrix is", signif(rcondAI , 2), "\n\tAI matrix may be singular - switching to an iteration of the EM algorithm\n")
      }


#TODO XXX XXX Need to return a message so can do em() outside ai(), but still within the current iteration (some exit message necessary, as part of class or structure) - maybe if AI all NAs on exit of the function then execute em()
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
   return(structure(list(thetav = thetavout, AI = AI, dLdtheta = dLdtheta),
	class = "gremlin"))

  }
################################################################################    











################################################################################
#' @rdname gremlinRmod
#' @export
#Meyer 1996:
## Likelihood eqn 3 <--> log determinants
## partial first & second derivatives of likelihood wrt parameters <--> eqn 4 and 5
##XXX This is broken down into sum of compoments: XXX
### log|R| <--> eqns 20, 25 & 26 (see univariate simplification in text following)
### log|G| <--> eqns 28, 31 & 32 (see univariate simplificaation in text following) 
### log|C| and tyPy: eqns 6, 9, 12-14 give how to get aspects of C and tyPy from M
#### log|C| <--> eqn 13
#### tyPy <--> eqn 14
#### first and second partial derivatives of log|C| and tyPy: eqn 15-18
#####XXX All above is before AI algorithm
aiNew <- function(thetavin, skel, thetaG, thetaR, sigma2e, modMats, W, sLc, sln){#Cinv, nminfrfx, r, ezero){
# thetavin <- thetav; sigma2e <- remlOut$sigma2e
    p <- length(thetavin)
#    thetain <- vech2matlist(thetavin, skel)
    AI <- matrix(NA, nrow = p, ncol = p)
    si <- modMats$nb+1
    B <- Matrix(0, nrow = modMats$ny, ncol = p, sparse = TRUE)
    for(g in thetaG){ #FIXME assumes thetaG is same length as thetavin
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # Meyer 1997 MS: eqn. 20 [also f(theta) in Johnson & Thompson 1995 eqn 9c)]
      # Knight 2008: eqn 2.32-2.35 and 3.11-3.13 (p.34)
      B[, g] <- modMats$Zg[[g]] %*% sln[si:ei, , drop = FALSE] %*% Ginv[[g]] #FIXME is order correct?
    }
#FIXME TODO Check what to do if more than 1 residual variance parameter
    if(g < p){
      B[, p] <- r / sigma2e
    }      
    iRHS <- Matrix(crossprod(W, B), sparse = TRUE)
    isln <- solve(a = sLc, b = iRHS, system = "A")
    # calculate residuals
#TODO are residuals y-W%*%isln OR bi-W%*%isln?????????????XXX
    PB <- B - W %*% isln
    AI <- (1 / 2*sigma2e) * as(crossprod(B, PB), "matrix")


#    if((g+1) < p){
#        for(k in (g+1):(p-1)){  #<-- fill upper triangle
#          sk <- sum(sapply(seq(k-1), FUN = function(u){ncol(modMats$Zg[[u]])})) + modMats$nb + 1
#          ek <- sk - 1 + ncol(modMats$Zg[[k]])
#          AI[g, k] <- AI[k, g] <- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% W[, sk:ek] %*% sln[sk:ek, , drop = FALSE]) / (thetavin[g]*thetavin[k]))@x    
#        }  #<-- end 'for(k ...)`
#      }  #<-- end `if()`
#      AI[g, p] <- AI[p, g] <- 0.5 * ((crossprod(isln, crossprod(W[, si:ei], PorVinv)) %*% r) / (thetavin[g]*thetavin[p]))@x
#      si <- ei+1
#    }   #<-- end `for(g ...)`
#    AI[p, p] <- 0.5 * ((crossprod(r, PorVinv) %*% r) / thetavin[p]^2)@x


return(list(AI = AI))
#   return(structure(list(thetav = thetavout, AI = AI, dLdtheta = dLdtheta),
#	class = "gremlin"))

  }
################################################################################    







################################################################################
#' @rdname gremlinRmod
#' @export
  gradFun <- function(thetav, thetaG, thetaR, modMats, Cinv, nminfrfx, sln, r){
    dLdtheta <- matrix(NA, nrow = p, ncol = 1, dimnames = list(names(thetav), NULL))
    # tee = e'e
    tee <- crossprod(r)
    # trCinvGeninv_gg = trace[Cinv_gg %*% geninv_gg] | tugug = u_gg' %*% geninv_gg %*% u_gg
    ## g is the gth component of the G-structure to model
    ## geninv is the generalized inverse (not the inverse of the G-matrix/varcomps)
#TODO make variable `length(thetaG)`
    trCinvGeninv_gg <- tugug <- as.list(rep(0, length(thetaG)))
    si <- modMats$nb+1
    for(g in thetaG){
      qi <- ncol(modMats$Zg[[g]])
      ei <- si - 1 + qi
      # note the trace of a product is equal to the sum of the element-by-element product
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      tugug[[g]] <- crossprod(sln[si:ei, , drop = FALSE], modMats$listGeninv[[g]]) %*% sln[si:ei, drop = FALSE]
#TODO FIXME check what to do if no ginverse associated with a parameter?!
      trCinvGeninv_gg[[g]] <- tr(modMats$listGeninv[[g]] %*% Cinv[si:ei, si:ei])
      si <- ei+1
    } 

    # First derivatives (gradient)
#TODO check that `nminfrfx` is still `n-p-q` when more than 1 random effect 
## ALSO check what `q` of `n-p-q` is when >1 random effect
#FIXME change `[p]` below to be number of residual (co)variances
    dLdtheta[p] <- 0.5*((tee / tail(thetav, 1)^2) - (nminfrfx)/tail(thetav, 1)) 
    for(g in thetaG){
      dLdtheta[p] <- dLdtheta[p] - 0.5*(trCinvGeninv_gg[[g]]/thetav[g]) 
      dLdtheta[g] <- 0.5*(tugug[[g]]/(thetav[g]^2) - ncol(modMats$Zg[[g]])/thetav[g] + trCinvGeninv_gg[[g]]*tail(thetav, 1)/(thetav[g]^2))
    }

   dLdtheta
  }
#############################
# Separating and pre-allocating P and Vinv to sparse Matrix doesn't seem to make
## much of a difference
# Also changing ginverse elements to `dsCMatrix` doesn't speedup traces, since
## they end up more or less as dense matrices but in dgCMatrix from the product
#############################

