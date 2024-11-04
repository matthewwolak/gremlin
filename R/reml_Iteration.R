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
#'   grSsetp <- gremlinSetup(WWG11 ~ sex - 1, random = ~ sire, data = Mrode11)
#'   grS <- remlIt(grSsetp)
#'
#' @export
###########
# Generic
###########
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
  ## Below also checks that all are implemented in cpp (if nomatch <--> NA)
  ### (currently only EM (1) or AItr/AIfd (2))
  intfacalgit <- rep(NA, grMod$maxit)
    intfacalgit[which(grMod$algit[1:grMod$maxit] == "EM")] <- 1
    intfacalgit[which(grMod$sdit[1:grMod$maxit] == "AI")] <- 2
    if(any(is.na(intfacalgit))){
      stop(cat("Algorithm", grMod$algit[which(is.na(intfacalgit))],
	"not implemented in c++, try", sQuote('gremlinR()'), "\n"))
    }
  intfacfdit <- as.integer(grMod$fdit)
  
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
      if(grMod$v > 0){  
        cat("\n\n\t\t\t **NOTE:**",
        "\ngremlin's c++ code is far too sophisticated for such a simple model.\n",
        "\tRe-fitting with the", sQuote('gremlinR()'), "function instead\n\n")
      }
      return(remlIt.gremlinR(grMod))
    }

  bounds <- grMod$bounds
    bounds[is.na(bounds)] <- 0  # Replace NAs with 0


  Cout <- .C("ugremlin", PACKAGE = "gremlin",
	as.double(grMod$modMats$y),
	as.integer(grMod$modMats$ny),
	as.integer(grMod$nminffx),		# No. observations - No. Fxd Fx
	as.integer(grMod$ndgeninv),		# non-diagonal ginverses
	as.integer(c(grMod$dimsZg)),
	as.integer(dimZWG), 			# Z, W, and geninv Dims
	as.integer(nnzWG), 			# No. nnz in W and geninvs
	as.integer(grMod$W@i), 			#W
	as.integer(grMod$W@p),
	as.double(grMod$W@x),
	as.integer(geninv_i), 			#geninv (generalized inverses)
	as.integer(geninv_p),
	as.double(geninv_x),
	as.double(grMod$rfxIncContrib2loglik),	# Random Fx contribution to log-Likelihood
        as.integer(grMod$lambda),	# TRUE/FALSE lambda (variance ratio?)
	as.integer(grMod$p),				#p=No. nu params
	as.integer(c(length(grMod$thetaG), length(grMod$thetaR))), #No. G and R nus
	as.integer(unlist(lapply(gnu, FUN = function(g) g@Dim[[1L]]))),#dim GRs
	as.integer(unlist(lapply(gnu, FUN = function(g) g@i))),	      #i GRs
	as.integer(unlist(lapply(gnu, FUN = function(g) g@p))),	      #p GRs
	as.integer(unlist(lapply(gnu, FUN = function(g) length(g@x)))), #no. non-zero GRs
	as.double(unlist(lapply(gnu, FUN = function(g) g@x))),	#nu vector
        as.integer(as.integer(grMod$conv)-1),   	# constraint codes (F=0)
        as.double(c(bounds)),			#boundaries (1:p=LB | p+1:2p=UB
	as.integer(length(grMod$Bpinv@x)),	#Bpinv (fixed fx prior inverse)
	as.integer(grMod$Bpinv@i),
	as.integer(grMod$Bpinv@p),
	as.double(grMod$Bpinv@x),
	as.double(rep(0, grMod$p)),			#empty dLdnu
	as.double(rep(0, grMod$p^2)),	#empty column-wise vector of AI matrix
	as.double(c(grMod$sln)), 			#empty sln
        as.double(c(grMod$Cinv_ii)),			#empty diag(Cinv)
	as.double(c(grMod$r)),				#empty resdiuals
	as.double(rep(0, grMod$maxit*(grMod$p+5))),	#itMat
	as.integer(intfacalgit - 1), 			#algorithm for each iteration
	as.integer(intfacfdit - 1),			#first deriv. algorithm
	as.integer(grMod$maxit),			#max it./n algit
	as.double(grMod$step),			#init./default step-halving value
	as.double(grMod$cctol),				#convergence tol.
	as.double(grMod$ezero),				#effective 0
        as.double(grMod$einf),				#effective +/- max. values
	as.integer(grMod$v),				#verbosity
	as.integer(grMod$vit),				#when to output status
	as.double(grMod$h),				#finite difference h value
	as.integer(rep(0, length(grMod$sln))))		#empty sLc->pinv

  i <- Cout[[37]]  #<-- index from c++ always increments +1 at end of for `i`

  grMod$nu[] <- vech2matlist(Cout[[22]], attr(grMod$thetav, "skel"))
  grMod$dLdnu[] <- Cout[[29]]
  if(all(Cout[[30]] == 0)) grMod$AI <- NULL else{
    grMod$AI <- matrix(Cout[[30]], nrow = grMod$p, ncol = grMod$p, byrow = FALSE)
      dimnames(grMod$AI) <- list(rownames(grMod$dLdnu), rownames(grMod$dLdnu))
  }

  grMod$conv[] <- levels(grMod$conv)[ Cout[[23]] + 1 ]
  grMod$sln[] <- Cout[[31]]
  grMod$Cinv_ii <- Cout[[32]] 
  grMod$r[] <- Cout[[33]]
  #TODO Will definitely need R vs. c++ methods for `update.gremlin()`
  #### can directly use R's `grMod$sLc`, but will need to figure out how to give c++'s `cs_schol()` a pinv (need to reconstruct `sLc` in c++ around pinv (see old code on how I may have done this when I made sLc from sLm)
  grMod$sLcPinv <- Cout[[45]]

  intfacalgitOut <- Cout[[35]][1:i] + 1
  if(any(intfacalgitOut > 2)) warning("unrecognized algorithm outputted")
  intfacfditOut <- Cout[[36]][1:i] + 1
  if(any(intfacfditOut > 5)){
    warning("unrecognized first derivative algorithm outputted")
  }
    intfacalgitEM <- which(intfacalgitOut == 1) #<-- index used several times
    intfacalgitAI <- which(intfacalgitOut == 2) #<-- index used several times
    # change second derivative values to NA for iterations AI fail, switch to EM
    grMod$sdit[1:i][intfacalgitEM] <- NA
    # change first derivative values for iterations AI fail, switch to EM
    grMod$fdit[1:i][intfacalgitEM] <- levels(grMod$fdit)[5]
    # remake algit values
    grMod$algit[1:i][intfacalgitEM] <- "EM"
    grMod$algit[1:i][intfacalgitAI] <- with(grMod,
      paste0(sdit[1:i][intfacalgitAI], fdit[1:i][intfacalgitAI]))
      
  itMat <- matrix(Cout[[34]][1:(i*(grMod$p+5))], nrow = i, ncol = grMod$p+5,
           byrow = TRUE)
    dimnames(itMat) <- list(paste(seq(i), grMod$algit[1:i], sep = "-"),
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
  conv <- grMod$conv
  bounds <- grMod$bounds
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
  prtCinv <- NULL
  prtCinvIt <- 0  #<-- keep latest iteration in which prtCinv was calculated

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
      sLc <- remlOut$sLc

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



    if(grMod$v > 0 && vitout == 0){
      cat("\tlL:", format(round(itMat[i, "loglik"], 6), nsmall = 6))
    }




    #################################
    # 5c check convergence criteria
    cc <- rep(NA, 4)
    if(i > 1){
      cc[1] <- ccFun1()   
      cc[2] <- ccFun2()
    } else cc[1] <- FALSE  #<-- ensures one of the EM/AI/etc algorithms used if i==1



    #################################
    # 5d Determine next (co)variance parameters to evaluate:
    ## if any algorithms need the diagonals of Cinv, calculate for this iteration
    if(grMod$fdit[i] == "tr" | grMod$algit[i] == "EM"){	
      if(grMod$v > 1 && vitout == 0) cat("\n\tcreate partial Cinv: chol2cinv_ii")
      prtCinv <- chol2inv_ii(expand(sLc)$L, Z = prtCinv)
        prtCinvIt <- i
      # Calculate trace and u_g'u_g (tugug) terms needed in EM and "tr"/gradFun
      ## only do this if no variance components besides residual
      if(length(thetaG) > 0){
        if(grMod$v > 1 && vitout == 0) cat("\n\tcalculate tugug(s) & trace(s)")
        tugug_tr_out <- tugug_trace(thetaG, grMod$modMats$nb, grMod$rfxlvls,
              grMod$modMats$listGeninv, prtCinv, grMod$sln,
              pinv = sLc@perm + as.integer(1)) 
      }
    }

    ############################
    #    EM
    ############################
    if(grMod$algit[i] == "EM" && !all(cc, na.rm = TRUE)){
      if(grMod$v > 1 && vitout == 0) cat("\n\tEM to find next nu")
      nuvout <- em(nuv, thetaG, thetaR, conv,
          grMod$rfxlvls, tugug_tr_out$tugug, tugug_tr_out$trace,
          grMod$modMats$y, grMod$r, grMod$nminffx)
    }


    ########################################
    #    AI with either analytical or
    ##  finite difference first derivatives
    ########################################
    if(grMod$algit[i] != "EM" && grMod$sdit[i] == "AI"){
      if(grMod$v > 1 && vitout == 0) cat("\n\tAI to find next nu")
#FIXME Currently, only allow when not: 
if(nrow(theta[[thetaR]]) != 1){
  stop(cat("\nAI algorithm currently only works for a single residual variance"))
}

      if(grMod$fdit[i] != "tr"){
        # finite difference algorithm for first derivatives
        #TODO instead interface next line with gremlinControl$algArgs
        if(grMod$fdit[i] == "cfd"){
          fd_in <- "cdiff"
        } else{
            fd_in <- ifelse(grMod$fdit[i] == "bfd", "bdiff", "fdiff")
          }  
        dLdnu <- gradFun_fd(nuvin = nuv, grObj = grMod,
          lL = remlOut$loglik, fd = fd_in)
      }  #<-- end finite difference first derivative choice

      if(lambda){
        AI <- ai(nuv, skel, thetaG,
             grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
  	      thetaR = NULL,
	      sigma2e)  #<-- NULL if lambda==FALSE
	if(grMod$fdit[i] == "tr"){  #<-- analytical first derivatives      
          dLdnu <- gradFun(nuv, thetaG, grMod$rfxlvls, grMod$sln,
          	      tugug_tr_out$tugug, tugug_tr_out$trace,
	    	      sigma2e = sigma2e, r = NULL, nminfrfx = NULL)
        } #<-- end analytical first derivative choice
          
      } else{
          AI <- ai(nuv, skel, thetaG,
        	grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
                thetaR,   #<-- NULL if lambda==TRUE
	        sigma2e = NULL)

	  if(grMod$fdit[i] == "tr"){  #<-- analytical first derivatives      
	    dLdnu <- gradFun(nuv, thetaG, grMod$rfxlvls, grMod$sln,
          	        tugug_tr_out$tugug, tugug_tr_out$trace,
  	      	        sigma2e = NULL, grMod$r, grMod$nminfrfx)

          } #<-- end analytical first derivative choice
        }

      ## Find next set of parameters using a quasi-Newton method/algorithm
      ### Meyer 1989 pp. 326-327 describes quasi-Newton methods 
#TODO see Meyer 1997 eqn 58 for Marquardt 1963: theta_t+1=theta_t - (H_t + k_t * I)^{-1} g_t 
## What I do below is similar: except k_t=f
      ### Mrode 2005 eqn 11.4
      ### Johnson and Thompson 1995 eqn 12
      ####(though gremlin uses `+` instead of J & T '95 `-` because
      ##### gremlin multiplies gradient by -0.5 in `gradFun()`)

      # Check for fixed (co)variance parameters
      ## remove them from copies of gradient and AI if so
      fxdP <- which(conv == "F")
      if(length(fxdP) > 0){
        AI_con <- AI[-fxdP, -fxdP, drop = FALSE]
        dLdnu_con <- dLdnu[-fxdP, , drop = FALSE]
      } else{
          AI_con <- AI
          dLdnu_con <- dLdnu
        }

      # Check if AI (after removing constrained parameters) can be inverted
      rcondAI <- rcond(AI_con)
      if(rcondAI < grMod$ezero){
        if(grMod$v > 2){
          cat("\n\tReciprocal condition number of H matrix is",
	    signif(rcondAI , 2), 
	    "\n\t   H matrix may be singular - modifying diagonals")
        }  #<-- end `if v>2`
      }  #<-- end if AI singular
      ### Check/modify AI matrix to 'ensure' positive definiteness
      ### `fI` is factor to adjust AI matrix
      #### (e.g., Meyer 1997 eqn 58 and WOMBAT manual A.5 strategy 3b)
      AIeigvals <- eigen(AI_con, symmetric = TRUE, only.values = TRUE)$values
        d <- (3e-6) * AIeigvals[1]
        f <- max(0, d - AIeigvals[nrow(AI_con)])
      fI <- f * diag(x = 1, nrow = nrow(AI_con))
      ##### modified 'Hessian'
      H <- fI + AI_con
      # Check if H can be inverted
      rcondH <- rcond(H)
      ## if H cannot be inverted do EM
      if(rcondH < grMod$ezero){
        if(grMod$v > 1){
          cat("\n\tReciprocal condition number of modified Hessian is",
	    signif(rcondH , 2), 
	    "\n\t   Hessian may be singular - switching to EM algorithm")
        }  #<-- end `if v>1`
        if(grMod$v > 1 && vitout == 0) cat("\n\tEM to find next nu")
          nuvout <- em(nuv, thetaG, thetaR, conv,
            grMod$rfxlvls, tugug_tr_out$tugug, tugug_tr_out$trace,
            grMod$modMats$y, grMod$r, grMod$nminffx)
          grMod$algit[i] <- "EM"
          grMod$fdit[i] <- "dfree"
          grMod$sdit[i] <- NA

      } else{  #<-- end if Hessian cannot be inverted
          Hinv <- solve(H)
          dnu <- Hinv %*% dLdnu_con   #<-- proposed change in nu parameters
          # First, implement step-reduction (if necessary)
 	  ## Rule: if `dnu` proposed greater than 200% change in any parameter 
          ### Then implement step reduction (`grMod$step` default) else do not
          nuvout <- matrix(nuv, ncol = 1)
          if(length(fxdP) > 0){
            if(any(abs(dnu / nuvout[-fxdP, ]) > 2.0)){
  	      step <- grMod$step
  	    } else step <- 1.0
            dnu <- step * dnu
            nuvout[-fxdP, ] <- nuvout[-fxdP, ] + dnu
            # Second, check for indecent proposals
            badLp <- sapply(seq(p)[-fxdP],
                  FUN = function(x) nuvout[x,] <= bounds[x, "LB"])
            badUp <- sapply(seq(p)[-fxdP],
                  FUN = function(x) nuvout[x,] >= bounds[x, "UB"])
            if(any(badLp) | any(badUp)){
              bad <- which((badUp + badLp) > 0)
              fxdPbadInd <- seq(p)[-fxdP][bad]
              if(grMod$v > 0){
                cat("\n(co)variance component(s) in 'thetav' vector (position:",
                    fxdPbadInd, ") restrained inside boundaries\t")
              }
              nuvout[-fxdP][which(badLp)] <- bounds[-fxdP, "LB"][which(badLp)] +
	          grMod$ezero
              nuvout[-fxdP][which(badUp)] <- bounds[-fxdP, "UB"][which(badUp)] -
                  grMod$ezero
              conv[fxdPbadInd] <- "B"  #<-- B for Bounded
            }  #<-- end if bad proposals AND fixed parameters

          } else{  #<-- now do below if NO fixed parameters
              if(any(abs(dnu / nuvout) > 2.0)){
                step <- grMod$step
              } else step <- 1.0
              dnu <- step * dnu
              nuvout <- nuvout + dnu
              # Second, check for indecent proposals
              badLp <- sapply(seq(p),
                    FUN = function(x) nuvout[x,] <= bounds[x, "LB"])
              badUp <- sapply(seq(p),
                    FUN = function(x) nuvout[x,] >= bounds[x, "UB"])
              if(any(badLp) | any(badUp)){
                bad <- which((badUp + badLp) > 0)
                if(grMod$v > 1){
                  cat("\n(co)variance component(s) in 'thetav' vector (position:",
                      bad, ") restrained inside boundaries\t")
                }
                nuvout[which(badLp), ] <- bounds[which(badLp), "LB"] +
	            grMod$ezero
                nuvout[which(badUp), ] <- bounds[which(badUp), "UB"] -
                    grMod$ezero
                conv[bad] <- "B"  #<-- B for Bounded
              }  #<-- end if bad proposals and NO fixed parameters
            }  #<-- end if/else fixed parameters


            ## Re-calculate parameter updates, CONDITIONAL on restrained values
            ### Gilmour 2019 AI REML in Practice. J. Anim. Breed. Genet.
            if(any(badLp) | any(badUp)){
              dLdnu_con[bad] <- 0.0  #<-- so convergence check 3 works correctly
              #TODO check in case all non-fixed parameters are bad!
              Hinv_uu <- solve(H[-bad, -bad, drop = FALSE])  #<-- un-restrained components
              H_uc <- H[, bad, drop = FALSE][-bad, , drop = FALSE]
              dLdnu_u <- dLdnu_con[-bad, , drop = FALSE]
              if(length(fxdP) > 0){
                nuvout[-c(fxdP, fxdPbadInd), ] <- nuv[-c(fxdP, fxdPbadInd)] +
                  (Hinv_uu %*% matrix((dLdnu_u - H_uc %*%
					(nuvout[fxdPbadInd]-nuv[fxdPbadInd])),
                    ncol = 1))
              } else{
                  nuvout[-bad, ] <- nuv[-bad] +
                    (Hinv_uu %*% matrix((dLdnu_u - H_uc %*% (nuvout[bad]-nuv[bad])),
                      ncol = 1))
                }
            } else{
                # Remove any boundary constraint codes from previous iterations
                ## restore to original code
                if(any(conv == "B")){
                  conv[which(conv == "B")] <- grMod$conv[which(conv == "B")]
                }
              }  #<-- end if/else indecent proposals

            # CONVERGENCE checks for AI
            cc[3] <- ccFun3()
#            cc[4] <- ccFun4()

        }  #<-- end if/else Hessian can be inverted


    }  #<-- end if algorithm is "AI"




    if(grMod$v > 1 && vitout == 0) cat("\n\tConvergence crit:", cc, "\n")





    ############################
    #    BOBYQA/NR
    ############################
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









    ############################################################################
    if(lambda) nuvout[thetaR] <- nuv[thetaR]  #<-- keep R=1 (R factored out)
    nuout <- vech2matlist(nuvout, skel) 
    nu <- sapply(nuout, FUN = stTrans)
    itTime <- Sys.time() - stItTime
    if(grMod$v > 0 && vitout == 0){
      if(grMod$v > 2){
        if(grMod$algit[i] != "EM"){
          sgd <- matrix(NA, nrow = p, ncol = p+2)  #<--`sgd` = summary.gremlinDeriv 
            dimnames(sgd) <- list(row.names(dLdnu),
              c("gradient", "", "AI", rep("", p-1)))
          sgd[, 1] <- dLdnu
#FIXME if get Hessian another way then need if statement about sdit          
          for(rc in 1:p) sgd[rc, 3:(p+2)] <- AI[rc, ]
          cat("\tstep reduction:", step, "\n")
          cat("\tH modification", round(f, 3), "\n")
          print(as.table(sgd), digits = 3, na.print = " | ", zero.print = ".")
          cat("\n")
        } 
      }  
      cat("\t\ttook ", round(itTime, 2), units(itTime), "\n")
    }
    units(itTime) <- "secs"
    itMat[i, ncol(itMat)] <- round(itTime, 1)



    # Determine if model has converged
    if(all(cc, na.rm = TRUE)){
      if(grMod$v > 0) cat("\n***  REML converged  ***\n\n")
      break
    }

  }  # END log-likelihood iterations
  #################################### 
  
  itMat <- itMat[1:i, , drop = FALSE]
    rownames(itMat) <- paste(seq(i), grMod$algit[1:i], sep = "-")
  if(lambda){
    theta <- nu2theta_lambda(nu, sigma2e, thetaG, thetaR)
  } else{
      theta <- nu2theta_noTrans(nu, thetaG, thetaR)
    }
  thetav <- matlist2vech(theta)


  # Calculate Cinv_ii and AI for last set of parameters
  ## Cinv_ii only if didn't calculate prtCinv for last iteration
  if(prtCinvIt < i){
    prtCinv <- chol2inv_ii(sLc, Z = prtCinv)
  }
  grMod$Cinv_ii[sLc@perm + as.integer(1)] <- prtCinv@x[attr(prtCinv, "Zdiagp")]
    
  ## AI
  if(grMod$sdit[i] != "AI" | is.na(grMod$sdit[i])){
    if(lambda){
      AI <- ai(nuv, skel, thetaG,
	     grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
	     thetaR = NULL,
	     sigma2e)  #<-- NULL if lambda==FALSE
    } else{ 
        AI <- ai(nuv, skel, thetaG,
   	     grMod$modMats, grMod$W, sLc, grMod$sln, grMod$r,
             thetaR,   #<-- NULL if lambda==TRUE
	     sigma2e = NULL)
      }
    dimnames(AI) <- list(rownames(dLdnu), rownames(dLdnu))
  }

  # place these altered values back into grMod
  grMod$conv <- conv
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



