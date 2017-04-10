# Matrix Trace methods
######################
###########
# Generic
tr <- function(X, ...){
  UseMethod("tr", X)
}

###########
# Default
tr.default <- function(X, ...){
  sum(diag(X))
}

###########
# dgCMatrix
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
tr.dsCMatrix <- function(X, ...){
  if(X@uplo == "L"){
    return(sum(X@x[(X@p[-(X@Dim[2L]+1)]+1)]))
  } else{
      return(sum(X@x[X@p[-1]]))
    }
}
  

