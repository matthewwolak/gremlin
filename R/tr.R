# Matrix Trace methods
######################
###########
# Generic
###########
#' Matrix trace methods.
#'
#' Methods to efficiently calculate a matrix trace depending on the class of matrix.
#'
#' @aliases tr tr.default tr.dgCMatrix tr.dsCMatrix
#' @param X A matrix.
#' @param \dots Additional arguments.
#'
#' @return A \code{numeric} value for the sum of the diagonal elements.
#' @author \email{matthewwolak@@gmail.com}
#' @export
tr <- function(X, ...){
  UseMethod("tr", X)
}

###########
# Default
###########
#' @describeIn tr Default method
#' @export
tr.default <- function(X, ...){
  sum(diag(X))
}

###########
# dgCMatrix
###########
#' @describeIn tr Method for matrix \code{X} of class Matrix:::dgCMatrix
#' @export
#' @import Matrix
tr.dgCMatrix <- function(X, ...){
  X <- forceSymmetric(X)
  if(X@uplo == "L"){
    return(sum(X@x[(X@p[-(X@Dim[2L]+1)]+1)]))
  } else{
      return(sum(X@x[X@p[-1]]))
    }
}

###########
# dsCMatrix
###########
#' @describeIn tr Method for matrix \code{X} of class Matrix:::dsCMatrix
#' @export
#' @import Matrix
tr.dsCMatrix <- function(X, ...){
  if(X@uplo == "L"){
    return(sum(X@x[(X@p[-(X@Dim[2L]+1)]+1)]))
  } else{
      return(sum(X@x[X@p[-1]]))
    }
}
  

