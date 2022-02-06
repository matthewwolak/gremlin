#include "gremlin.h"

// invert square matrix, possibly with fixed/constraints
//// remove constrained portions then check if can be inverted first
////// if issues with diagonals then add small amount and try again

/* inverted matrix returned if successful else NULL */
cs *cs_inv_withDiagMod(const cs *A, csi *con, csi *wchBd,
	double *ezero, csi v
){

  cs	*H, *Hinv;
  css	*sLh;
  csn	*Lh;
  
  csi	g, k, p, conP;
  
  double d, f = 0.0;
  
  if(!CS_CSC(A) ) return (NULL);    // check arguments

  conP = p = A->n;

  // Check for fixed (co)variance parameters
  for(k = 0; k < p; k++){
    if(con[k] == 0 || con[k] == 3){
      wchBd[k] = 1;  // only reset elements about to use
      conP--;
    } else wchBd[k] = 0;
  }

  H = cs_droprowcol(A, wchBd);      

  // Check if H can be inverted
  sLh = cs_schol(1, H);
  Lh = cs_chol(H, sLh);
  if(Lh == NULL){
    f = 3e-5; 

    if(v > 1){
      Rprintf("\n\tH cholesky decomposition failed:\n\t   Hessian matrix may be singular - modifying diagonals and re-trying");
        if(v > 2) Rprintf("\n\tH modification: %6.3g\n", f);
    } // end if v>1
        
  }  // end if Lh=NULL


  // Check/modify H matrix to 'ensure' positive definiteness
  /* `f` is factor to adjust H matrix
     (e.g., based off Meyer 1997 eqn 58 and WOMBAT manual A.5 strategy 3b)
     HOWEVER, differs from above in that I do not use eigenvalues here,
     rather the diagonals of the Cholesky factorization are inspected
       see `remlIt.gremlinR()` for approach like Meyer (with eigenvalues) */
/*
See `La_rs` defined around line 153 of "R/src/modules/lapack/Lapack.c"
What R's `eigen()` calls
 */

  //// check if any small/zero diagonals on Cholesky of H (Lh)
  ////// If any diagonals of Lh < 0 modify H by adding f to diagonals
  if(Lh != NULL){
    for(g = 0; g < Lh->L->n; g++){
      d = Lh->L->x[ Lh->L->p[g] ];
      if(d < ezero[0]){
        f = 3e-5;
        if(v > 2) Rprintf("\n\tSmall diagonal on Cholesky of H: adding %6.3g\n", f);

      }
    }
  }
  ////XXX ASSUME H is full matrix so H->p[g] + g = diagonal
  for(g = 0; g < H->n; g++) H->x[ H->p[g] + g ] += f;

  // Now check AGAIN if Hessian can be inverted
  cs_sfree(sLh); 
  cs_nfree(Lh);
  sLh = cs_schol(1, H);
  Lh = cs_chol(H, sLh);
  if(Lh == NULL){
    if(v > 1){
      Rprintf("\n\tH cholesky decomposition failed:\n\t   Hessian matrix may be singular");
    } // end if v>1

    cs_sfree(sLh);
    cs_nfree(Lh);
   return (cs_spfree(H));

  } //<-- end if H cannot be inverted  


  Hinv = cs_inv(H);

  // Cleanup:
  cs_sfree(sLh); 
  cs_nfree(Lh);
  cs_spfree(H);
 
 return (Hinv) ;	/* success; return Hinv */
}
  
