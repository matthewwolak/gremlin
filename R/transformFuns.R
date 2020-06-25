################################################################################
#' Transformation of starting parameters.
#'
#' Transform start parameters into format \code{gremlin} expects.
#'
#' @aliases stTrans
#' @param x A \code{list} of starting parameters.
#' @return A sparse \sQuote{dsCMatrix}
#' @author \email{matthewwolak@@gmail.com}
#' @export
#' @import Matrix
stTrans <- function(x){
  if(is.numeric(x) && !is.matrix(x)) x <- as.matrix(x)
  xnonames <- unname(x)
  if(!isSymmetric(xnonames)) stop(cat("Element", x, "must be a symmetric matrix or a number\n")) 
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
#' @export
#' @import Matrix
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



