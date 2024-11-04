#' Partial sparse matrix inverse from a Cholesky factorization.
#'
#' Only calculate values of a sparse matrix inverse corresponding to non-zero
#'   locations for the Cholesky factorization. 
#'
#' If $\strong{L L'} = \strong{C}$, function efficiently gives diag(Cinv) by only
#'   calculating elements of Cinv based on non-zero elements of $\strong{L}$ and
#'   $\strong{L}$. Follows the method and equations by Takahashi et al. (1973).
#'
#' @aliases chol2inv_ii
#' @param L A lower-triangle Cholesky factorization ($\strong{L L'} = \strong{C}$).
#' @param Z A sparse matrix containing the partial inverse of $\strong{L L'}$
#'   from a previous call to the function. Must contain the \dQuote{Zdiagp}
#'   attribute.
#'
#' @return A sparse matrix containing the partial inverse of C ($\strong{L L'}$)
#'   along with attribute \dQuote{Zdiagp} indicating the location for diagonals
#'   of Z in \code{slot x} of the object \code{Z}.
#'
#' @author \email{matthewwolak@@gmail.com}
#' @references
#'   Takahashi, Fagan, & Chin. 1973. Formation of a sparse bus impedance matrix
#'   and its application to short circuit study. 8th PICA Conference Proceedings,
#'   Minneapolis, MN.  
#' @export
#' @import Matrix
chol2inv_ii <- function(L, Z = NULL){
  n <- L@Dim[1L]
  Lp <- L@p + 1
  Li <- L@i + 1
  Lx <- L@x
  
  d <- 1 / Lx[ Lp[1:n] ]  #<-- diagonal of L
  
  ###############################
  if(is.null(Z)){
    # Find non-zero pattern of Z
    ## Zpattern = non-zero pattern of L + L'
    ## Figure this out without explicitly making L' or adding to L

    # calculate Zp
    ## start with L and tL column counts
    ### make tL WITHOUT diagonal
    tL <- triu(t(sparseMatrix(i = Li-as.integer(1), p = Lp-as.integer(1),
                 index1 = FALSE, symmetric = FALSE)), 1)
    tLi <- tL@i + as.integer(1)
    tLp <- tL@p + as.integer(1)
    Zp <- cumsum(c(1, diff(tLp) + diff(Lp)))
   
    # Add row values for Zi AND Find diagonal x of Z AND initialize Z[j, j] value
    ## initialize Zi, Zp, and Zx
    nnz <- Zp[n + 1] - 1
    Zx <- double(0)
    Zi <- integer(nnz)
    Zdiagp <- integer(n)
    
    pZp <- 1
    for(j in 1:n){
      # fill in upper triangle rows (not including diagonal)
      if(tLp[j] < tLp[j + 1]){
        pseq <- seq.int(tLp[j], tLp[j + 1] - 1, by = 1)  
        pZpseq <- seq.int(from = pZp, by = 1, along.with = pseq)        
          Zi[pZpseq] <- tLi[pseq] 
          pZp <- pZpseq[length(pseq)] + as.integer(1)
      }  #<-- end if
           
      pdiag <- -1  #should be >0; this will ID Z_jj/diagonal as zero=not allowed
      # Go through each row of L[, j] which is diagonal and below rows of Z[, j]
      pseq <- seq.int(Lp[j], Lp[j+1] - 1, by = 1)
      pZpseq <- seq.int(from = pZp, by = 1, along.with = pseq)
        # find when row index = j
        if(any(Lijdiag <- Li[pseq] == j)){
          pdiag <- pZpseq[Lijdiag]
        }
      Zi[pZpseq] <- Li[pseq] 
      pZp <- pZpseq[length(pseq)] + as.integer(1)
     
      if(pdiag == -1) return(-1) #<-- Z must have zero-free diagonal
      Zdiagp[j] <- pdiag 

    }  #<-- end for j
    
  } else{  #<-- end if Z=NULL
      Zi <- Z@i + as.integer(1) 
      Zp <- Z@p + as.integer(1)
      Zdiagp <- attr(Z, "Zdiagp") #<-- location of diagonals in Zi
      if(is.null(Zdiagp)){
        stop("'Z' in chol2inv_ii() must have non-NULL 'Zdiagp' attribute")
      }
    }
  ###############################

    
  Zx <- double(Zp[n+1]-1)     
  Zx[Zdiagp] <- d * d  #<--using LL' (if LDL', then = d[j] or maybe 1/d)
  z <- double(n)  #<-- zero on out - holds Z[, j] workspace
  Lmunch <- Lp[2:(n+1)]-1  #<-- points to last (working) entry in column k of L


  # Compute sparse inverse subset
  for(j in seq.int(n, 1, -1)){

    # scatter Z[, j] into z workspace
    ## only lower triangular part is needed, since upper is all zero to start
    ###TODO if only hold lower part Z this could be simpler (no need for Zdiagp)
    pseq <- seq(Zdiagp[j], Zp[j + 1] - 1, 1)
      z[ Zi[pseq] ] <- Zx[pseq]



    #XXX for R needs if check below
    if((Zdiagp[j] - 1) >= Zp[j]){
    for(p in seq(Zdiagp[j] - 1, Zp[j], by = -1)){
      # compute Z[, j] **above** diagonals
      ## for k=(j-1) to -1 by 1, but only for entries Z[k, j]
      ### moves from just above diagonal up to top row of Z[, j]

      # Z[k, j] = - U[k, (k+1):n] * Z[(k+1):n, j]
      k <- Zi[p]  #<-- upper triangle row of Z
      
      # U=t(L) and if U stored by row this is same as L stored by column (so L=U)
      ## at kth row Z[, j] starts in kth row of U (col of L) moves to end of U(L)
      zkj <- function(up){
        # skip diagonal of U
        i <- Li[up]
        if(i > k){
          # need to adjust L (LL') so it is L from LDL' by * (1/D[k])
          return( -1 * Lx[up] * d[k] * z[i])
        } else return(0)
                
      }  #<-- end zkj function for up
      # zkj <- zkj - (Lx[up] * d[k] * z[i])
      z[k] <- sum(vapply(seq(Lp[k], Lp[k+1]-1, by = 1), FUN = zkj, double(1)))
      

      # Left-looking update to Z **below** diagonals 
      ## for k=(j-1) to -1 by 1, but only for entries Z[k, j]
      ### moves from just above diagonal up to top row of Z[, j]
      # ljk = L[j, k]
      if(Lmunch[k] < Lp[k] | Li[ Lmunch[k] ] != j){
        next #<-- skip L[j, k] = 0 since no math/work to do
      }
      
      # need to adjust L (LL') so it is L from LDL' by * (1/D[k])
      ljk <- Lx[ Lmunch[k] ] * d[k]
      Lmunch[k] <- Lmunch[k] - as.integer(1)
      
      # Z[(k+1):n, k] = Z[(k+1):n, k] - Z[(k+1):n, j] * L[j, k]
      ## moves from diagonal of Z[, k] to last row of Z[, k]
      zpseq <- seq.int(Zdiagp[k], Zp[k+1] - as.integer(1), by = 1)
      Zx[zpseq] <- Zx[zpseq] - z[ Zi[zpseq] ] * ljk
      
    }  #<-- end for p
    }  #<-- end if statement XXX just for R
    
  
    
    # Gather Z[, j] back from z workspace
    ## moves down column Z[, j] by (non-zero) rows
    pseq <- seq.int(Zp[j], Zp[j + 1] - as.integer(1), by = 1)
      Zx[pseq] <- z[ Zi[pseq] ]
      z[ Zi[pseq] ] <- 0.0
      
  }  #<-- end for j
  
  
  return(structure(sparseMatrix(i = Zi-1, p = Zp-1, x = Zx,
             index1 = FALSE, symmetric = FALSE),
           Zdiagp = Zdiagp))

}  #<-- end chol2inv_ii





