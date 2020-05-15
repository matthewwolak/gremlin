#include "gremlin.h"


/* matrix inversion from Cholesky Decomposition: L t(L)=A and t(Linv) Linv = Ainv */
/* iterate through columns j:
forward solve with identity matrix: L b_j = i_j (i_j is column of identity matrix)
backwards solve:  Lt Ainv_j = b_j (Lt is transpose of L)
*/
/* replaces Cinv_ii with diagonals from inverted matrix  */
// `r` tells what row to start from so can only get diagonals of part of matrix

csi cs_chol2inv_ii(const cs *L, const csi *Pinv, double *Cinv_ii, int r){  

  csi	 i, j, k, n, *Lp, *Li;
  double *b, *Lx;

  n = L->n;
  Lp = L->p; Li = L->i; Lx = L->x;
  b = cs_malloc(n, sizeof (double));

  for(k = r; k < n; k++){
    for(i = 0; i < n; i++) b[i] = 0.0;  // clear b
    // initiate jth column of identity matrix  
    b[Pinv[k]] += 1.0;           /* essentially `cs_ipvec` */
    // forward solve (e.g., cs_lsolve) b = L\b
    for(j = Pinv[k]; j < n; j++){
      if(b[j] != 0.0){
        b[j] /= Lx[Lp[j]];  // set diagonal ( 1 / L[k,k])
        // for loop to determine off-diagonal contributions
        for(i = Lp[j]+1; i < Lp[j+1]; i++){
          b[Li[i]] -= Lx[i] * b[j];
        }
      }  // end if NOT 0
    }
 
    // back solve (e.g., cs_ltsolve) b = L' \ b
    //// only need to find diagonal element for Cinv_ii[k] so stop at diagonal
    for(j = n-1; j >= Pinv[k]; j--){
      for(i = Lp[j]+1; i < Lp[j+1]; i++){
        b[j] -= Lx[i] * b[Li[i]];
      }
      b[j] /= Lx[Lp[j]];
    }

    // load diagnoal from b into output
    Cinv_ii[k] = b[Pinv[k]];    /* essentially `cs_pvec` */

  }  // end for k
  cs_free(b);

 return(1);
}


