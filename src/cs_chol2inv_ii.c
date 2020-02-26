#include "gremlin.h"


/* matrix inversion from Cholesky Decomposition: L t(L)=A and t(Linv) Linv = Ainv */
/* iterate through columns j:
forward solve with identity matrix: L b_j = i_j (i_j is column of identity matrix)
backwards solve:  Lt Ainv_j = b_j (Lt is transpose of L)
*/
/* replaces Cinv_ii with diagonals from inverted matrix  */
// `r` tells what row to start from so can only get diagonals of part of matrix

csi cs_chol2inv_ii(const cs *L, const csi *Pinv, double *Cinv_ii, int r){  

  csi	 g, j, n;
  double *b;

  n = L->n;
  b = cs_malloc(n, sizeof (double));
  for(g = 0; g < n; g++) b[g] = 0.0;

  for(j = r; j < n; j++){
    // initiate jth column of identity matrix  
    b[Pinv[j]] = 1.0;
    cs_lsolve(L, b);    	/* b = L\b   */
    cs_ltsolve(L, b);   	/* b = L'\b  */
    cs_pvec(Pinv, b, b, n);     /* b = P'*x  */

    // load diagnoal from b into output
    Cinv_ii[j] = b[j];
Rprintf("\nb[%i]=%6.4f and b[%i]=%6.4f", j, b[j], j+1, b[j+1]);
    // Also, reset b for next column
    for(g = 0; g < n; g++) b[g] = 0.0;
  }  // end for j
  cs_free(b);

 return(1);
}


