# Knight 2008 (ch. 6) says Searle et al. 1992 and Longford 1993 discuss
## different types of converg. crit.
## See Appendix 2 of WOMBAT help manual for 4 convergence criteria used
#
#' Convergence Criteria Checks for REML.
#'
#' Determine whether the optimization has converged on a maximum of the 
#' log-likelihood function
#'
#' @aliases ccFun, ccFun1, ccFun2, ccFun3, ccFun4
#' @param obj Optional gremlin model object. If \code{NULL} then the necessary
#'   variables are taken from the parent environment, if present
#'
#' @return A \code{logical} value whether the current REML iteration has passed
#'   the convergence criteria
#' @references
#'   Meyer, K. 2019. WOMBAT A program for mixed model analyses by restricted
#'   maximum likelihood. User Notes. 27 September 2019.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#'   grS2 <- gremlinR(WWG11 ~ sex, random = ~ sire, data = Mrode11, maxit = 2)
#'     ccFun1(grS2)
#'     ccFun2(grS2)
#'   grS <- gremlinR(WWG11 ~ sex, random = ~ sire, data = Mrode11)
#'     ccFun(grS)
#'
#' @export
ccFun <- function(obj = NULL){
  if(is.null(obj)){
    convCrit <- c(ccFun1(), ccFun2(), ccFun3())#, ccFun4())
  } else convCrit <- c(ccFun1(obj), ccFun2(obj), ccFun3(obj))#, ccFun4(obj))

 if(all(convCrit)) return(TRUE) else return(FALSE)
}

################################################################
# wombat 1
#' @rdname ccFun
#' @export
ccFun1 <- function(obj = NULL){
  if(is.null(obj)){
    itMat <- get("itMat", parent.frame())
    cctol <- get("grMod", parent.frame())$cctol
    i <- get("i", parent.frame())
  } else{
      itMat <- obj$itMat
      cctol <- obj$grMod$cctol
      i <- nrow(itMat)
    }      
  
 as.vector(diff(itMat[c(i-1, i), "loglik"])) < cctol[1]
}


################################################################
# wombat 2 (eqn. A.1) (also Knight 2008 (eqn. 6.1) criteria
#' @rdname ccFun
#' @export
ccFun2 <- function(obj = NULL){
  if(is.null(obj)){
    itMat <- get("itMat", parent.frame())
    cctol <- get("grMod", parent.frame())$cctol
    i <- get("i", parent.frame())
    p <- get("grMod", parent.frame())$p
  } else{
      itMat <- obj$itMat
      cctol <- obj$grMod$cctol
      i <- nrow(itMat)
      p <- obj$grMod$p
    }      
  
 sqrt(sum((itMat[i, 1:p] - itMat[(i-1), 1:p])^2)/sum(itMat[i, 1:p]^2)) < cctol[2]
}

################################################################
# ccFun3 and ccFun4 are for AI algorithms only
## wombat 3 (eqn. A.2): Norm of the gradient vector
#' @rdname ccFun
#' @export
ccFun3 <- function(obj = NULL){
  if(is.null(obj)){
    itMat <- get("itMat", parent.frame())
    cctol <- get("grMod", parent.frame())$cctol
    i <- get("i", parent.frame())
    dLdnu <- get("grMod", parent.frame())$dLdnu
    conv <- get("grMod", parent.frame())$conv
    fdit <- get("grMod", parent.frame())$fdit
  } else{
      itMat <- obj$itMat
      cctol <- obj$grMod$cctol
      i <- nrow(itMat)
      dLdnu <- obj$grMod$dLdnu
      conv <- obj$grMod$dLdnu
      fdit <- obj$grMod$fdit
    }      
  # do checks if first derivative algorithm used, otherwise can't/don't
  if(is.na(fdit[i])){
    return(TRUE)
  } else{
      grad <- dLdnu[which(conv == "F")]
      return(sqrt(sum(grad * grad)) < cctol[3])
    }
}

################################################################
# wombat 4 (eqn A.3): Newton decrement
## (Boyd & Vandenberghe 2004 "Convex Optimization" book cited in wombat)
# AI only
#TODO: what is criteria? Not stated in Wombat see Boyd & Vandenberghe 2004
#' @rdname ccFun
#' @export
ccFun4 <- function(obj = NULL){
  if(is.null(obj)){
    cctol <- get("grMod", parent.frame())$cctol
    #TODO: figure it whether to use "_con" versions or not
    dLdnu <- get("dLdnu_con", parent.frame())
    invH <- get("H_con", parent.frame())
  } else{
      cctol <- obj$grMod$cctol
      dLdnu <- obj$grMod$dLdnu
      invH <- solve(obj$grMod$AI)
    }          
 return((-1 * c(crossprod(dLdnu, invH) %*% dLdnu)) < cctol[4])
}


