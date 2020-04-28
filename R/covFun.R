################################################################################
#' (Co)variance parameter transformations.
#'
#' Converts lists of (co)variance parameters either between \code{list} and 
#'   \code{vector} format or between the theta and nu scales.
#'
#' \itemize{
#'   \item{stTrans }{Transform start parameters into lower triangle matrices of
#'     class \code{dsCMatrix}.}
#'   \item{start2theta }{Converts lists of starting values for (co)variance
#'     parameters to a theta object used to structure the (co)variance components
#'     within gremlin.}
#'   \item{matlist2vech }{Converts a \code{list} of (co)variance parameter
#'     matrices to a vector with a \dQuote{skel} attribute.}
#'   \item{vech2matlist }{Converts a vector of (co)variance parameters to a list of
#'     covariance matrices.}
#'   \item{theta2nu_trans }{Transforms theta to nu scale by taking the Cholesky
#'     factor of each covariance matrix and then replacing the diagonals with 
#'     their (natural) logarithms. Done to ensure matrices are positive definite.}
#'   \item{nu2theta_trans }{Back transformation from \code{theta2nu_trans}:
#'     exponentiates the diagonal elements of each matrix then calculates the
#'     cross-product.}
#'   \item{theta2nu_lambda }{Transformation that factors out a residual variance
#'     so that \code{nu} contains the \sQuote{lambda} parameterization: ratios
#'     of variance parameters with the residual variance.}
#'   \item{nu2theta_lambda }{Back transformation from \code{theta2nu_lambda}.}
#'   \item{nu2theta_noTrans }{Structures \code{theta} when not transformed.}
#' }
#'
#' @aliases stTrans vech2matlist start2theta matlist2vech theta2nu_trans
#'   nu2theta_trans theta2nu_lambda nu2theta_lambda nu2theta_noTrans
#' @param x,theta,nu A \code{list} of matrices containing the (co)variance
#'       parameters of the model.
#' @param vech A \code{vector} of (co)variance parameters.
#' @param skeleton An example structure to map \code{vech} onto.
#' @param Gstart,Rstart A \code{list} of starting (co)variance values for the
#'   G-structure (random effects terms) or R-structure (residual).
#' @param name An (optional) character \code{vector} containing the (co)variance
#'   component names.
#' @param thetaG,thetaR A \code{vector} indexing the G-structure or R-structure
#'   components, respectively.
#' @param sigma2e A \code{numeric} estimate of the factored out residual
#'   variance from the mixed model equations (i.e., the \sQuote{lambda} scale)
#'   \eqn{\sigma^{2}_{e}}.
#' @return Functions are specified to mostly return either a \code{list} of
#'   matrices (structure as defined by the \dQuote{skel} attribute or in
#'   the \code{skeleton} object) or a \code{vector} containing the (co)variance
#'   parameters of the model. Additional list elements returned can be:
#'   \describe{
#'     \item{thetaG }{A \code{vector} indexing the G-structure components.}
#'     \item{thetaR }{A \code{vector} indexing the R-structure components.}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @import Matrix
#' @examples
#'   # User-specified starting parameters
#'   thetaOut <- start2theta(Gstart = list(matrix(1), matrix(2)),
#'     Rstart = matrix(3))
#'   ## convert to a vector and then back into a matrix list
#'   thetav <- matlist2vech(thetaOut$theta)
#'   theta <- vech2matlist(thetav, attr(thetav, "skel"))
#'     identical(thetaOut$theta, theta)  #<-- should be TRUE
#'   # lambda parameterization transformation
#'   nu <- theta2nu_lambda(theta, thetaOut$thetaG, thetaOut$thetaR)
#'   # back-transform from (lambda scale) nu to theta
#'   ## For example, when the sigma2e estimate=0.5
#'   theta2 <- nu2theta_lambda(nu, sigma2e = 0.5, thetaOut$thetaG, thetaOut$thetaR)
#' @name covFun
NULL




# Transformation of starting parameters to correct matrix format.
#' @rdname covFun
#' @export
stTrans <- function(x){
  if(is.numeric(x) && !is.matrix(x)) x <- as.matrix(x)
  if(!isSymmetric(x)) stop(cat("Element", x, "must be a symmetric matrix or a number\n")) 
  x <- as(x, "symmetricMatrix")
  x@uplo <- "L"
  x <- as(x, "dsCMatrix")
 x
}



# Starting Parameters to theta List.
#' @rdname covFun
#' @export
start2theta <- function(Gstart, Rstart, name = NULL){
  theta <- c(G = sapply(Gstart, FUN = stTrans), R. = stTrans(Rstart))
  thetaGorR <- sapply(strsplit(names(theta), ".", fixed = TRUE), FUN = "[[", i = 1)
#FIXME ensure grep() is best way and won't mess up when multiple Gs and **Rs**
  thetaG <- grep("G", thetaGorR, ignore.case = FALSE)
  thetaR <- grep("R", thetaGorR, ignore.case = FALSE)
    thNames <- if(!is.null(name)) name else seq_len(length(thetaG))
    names(theta) <- c(paste0(rep("G.", length(thetaG)), thNames),
                      paste0("ResVar", seq(length(thetaR))))

 return(list(theta = theta, thetaG = thetaG, thetaR = thetaR))
}





# theta List to Vector
#' @rdname covFun
#' @export
matlist2vech <- function(theta){
  thetav <- sapply(theta, FUN = slot, name = "x")
  skel <- lapply(seq(length(theta)),
	FUN = function(i){mapply(slot, theta[i], c("i", "p", "Dim"))})
    names(skel) <- names(theta)
  attr(thetav, "skel") <- skel

 thetav
}


# Vector to List of Matrices.
#' @rdname covFun
#' @export
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







############################
# Log-diagonals of Cholesky
############################
#FIXME: change G to cholesky of G with log(diagonals)
### when to do cholesky factorization of G?
#### is it easier to do direct product of cholesky's than direct product then cholesky?
#### If the latter then save the symbolic factorization and just do updates?
#XXX Don't need for EM algorithm
#TODO FIXME need to take log of just diagonals
#TODO also convert to lambda scale

# Transform theta to nu With log-Cholesky Parameterization
#' @rdname covFun
#' @export
#' @import Matrix
theta2nu_trans <- function(theta){
  list(G = lapply(theta$G, FUN = function(x){log(chol(x))}),
       R = log(chol(theta$R)))
}

# Back-transform nu to theta, Reversing log-Cholesky Parameterization
#' @rdname covFun
#' @export
#' @import Matrix
#TODO FIXME when converting back, do exp() of just diagonals
#TODO also back-convert from lambda scale
nu2theta_trans <- function(nu){
  list(G = lapply(nu[[1L]], FUN = function(x){exp(crossprod(x))}),
       R = exp(crossprod(nu[[2L]])))
}



######################
# Lambda scale
######################

# Transform theta to nu on the lambda Scale
#' @rdname covFun
#' @export
#' @import Matrix
theta2nu_lambda <- function(theta, thetaG, thetaR){
  nu <- theta
  cRinv <- solve(chol(theta[[thetaR]]))
  # Meyer 1991, p.77 (<-- ?p70/eqn4?)
  nu[thetaG] <- lapply(thetaG,
    FUN = function(x){as(as(crossprod(cRinv, nu[[x]]) %*% cRinv, "symmetricMatrix"),
      "dsCMatrix")}) 
  #TODO what to do if R is a matrix
  nu[[thetaR]] <- as(as(crossprod(cRinv, nu[[thetaR]]) %*% cRinv, "symmetricMatrix"),
    "dsCMatrix")

 nu
}

# Back-transform lambda Scale nu to theta
#' @rdname covFun
#' @export
#' @import Matrix
nu2theta_lambda <- function(nu, sigma2e, thetaG, thetaR){ #TODO FIXME
  cR <- chol(matrix(sigma2e))
  theta <- nu
  theta[thetaG] <- lapply(thetaG,
    FUN = function(x){as(as(crossprod(cR, nu[[x]]) %*% cR, "symmetricMatrix"),
      "dsCMatrix")})
  theta[[thetaR]] <- as(as(crossprod(cR, nu[[thetaR]]) %*% cR, "symmetricMatrix"),
    "dsCMatrix")

 theta
}  



######################
# No transformations
######################

## Mostly ensures correct formatting of theta
## TODO combine with above function (reconcile with `vech2matlist()`)
#' @rdname covFun
#' @export
#' @import Matrix
# theta2nu_noTrans <- function(theta) theta

# Back-transform nu to theta When No Transformation Occurred
nu2theta_noTrans <- function(nu, thetaG, thetaR){
  theta <- nu
  theta[thetaG] <- lapply(thetaG,
    FUN = function(x){as(as(matrix(nu[[x]],1,1), "symmetricMatrix"),
      "dsCMatrix")})
  theta[[thetaR]] <- as(as(matrix(nu[[thetaR]], 1, 1), "symmetricMatrix"),
    "dsCMatrix")

 theta
}


