#include "gremlin.h"


/* matrix inversion from Cholesky Decomposition: L t(L)=A and t(Linv) Linv = Ainv */
/* iterate through columns j:
        forward solve with identity matrix: L b_j = i_j (i_j is column of identity matrix)
        backwards solve:  Lt Ainv_j = b_j (Lt is transpose of L)
*/
/* Gives lower triangle of Ainv so store as triplet matrix and then compress  */


cs *cs_chol2inv(const cs *L){  

  csi	 g, j, cnt, n, *Ainvi, *Ainvj;
  cs	 *Ainv, *INV;
  double *Ainvx, *b;

  n = L->n;
  b = cs_malloc(n, sizeof (double));
  Ainv = cs_spalloc(n, n, n*n, 1, 1);
  if(!CS_CSC (L) || !Ainv){
    cs_spfree(Ainv);
    return(NULL);
  } 
  Ainvi = Ainv->i; Ainvj = Ainv->p; Ainvx = Ainv->x;

  for(g = 0; g < n; g++) b[g] = 0.0;
  cnt = 0;
  //TODO if only want part of Inverse matrix, adjust which columns calculated
  for(j = 0; j < n; j++){
    // initiate jth column of identity matrix  
    b[j] = 1.0;
    cs_lsolve(L, b);    /* b = L\b  */
    cs_ltsolve(L, b);   /* b = L'\b */

    // load elements of b as Ainv[, j] and Ainv[j, ]
    //// Also, reset b for next column
    for(g = j; g < n; g++){
      Ainvi[cnt] = g; Ainvj[cnt] = j; Ainvx[cnt] = b[g];
      cnt++;
      if(g > j){
        Ainvi[cnt] = j; Ainvj[cnt] = g; Ainvx[cnt] = b[g];
        cnt++;
      }
      b[g] = 0.0;
    }  // end for g
    // finish resetting b
    if(j > 0) for(g = 0; g < j; g++) b[g] = 0.0;
  }  // end for j
  Ainv->nz = cnt - 1;

  INV = cs_compress(Ainv);    
  cs_spfree(Ainv);
  cs_free(b);
 return(INV);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



/* matrix inversion from Cholesky Decomposition: L t(L)=A and t(Linv) Linv = Ainv */
// See the paper 'Matrix Inversion Using Cholesky Decomposition'
//// by Aravindh Krishnamoorthy & Deepak Menon
//// arXiv:1111.4144.
/*

cs *cs_chol2inv(const cs *L){  

  csi	 g, i, j, k, *Li, *Lp, *Ainvi, *Ainvp, n, n2;
  cs	 *Ainv;
  double *Lx, *Ainvx, Ainvij, Sij, raij;

  n = L->n;
  n2 = n*n;
  Lp = L->p; Li = L->i; Lx = L->x;
  Ainv = cs_spalloc(n, n, n2, 1, 0);
  Ainv->p[0] = 0;
  for(k = 0; k < n; k++){
    Ainv->p[k+1] = Ainv->p[k] + n;
    for(i = 0; i < n; i++){
      Ainv->i[Ainv->p[k]+i] = i;
      Ainv->x[Ainv->p[k]+i] = 0.0;
    }
  }
  if(!CS_CSC (L) || !Ainv){
    cs_spfree(Ainv);
    return(NULL);
  } 

  Ainvp = Ainv->p; Ainvi = Ainv->i; Ainvx = Ainv->x;
  for(j = n-1; j >= 0; j--){
Rprintf("j=%i\n", j);
    for(i = j; i >= 0; i--){
Rprintf("\ti=%i\n", i);
      if(i == j) Sij = 1.0 / Lx[Lp[j]]; else Sij == 0.0;
      raij = 0.0;
      // raij results in answer to a [1 x n] * [n x 1] matrix multiplication
      // multiply each row in L[, i] to the sub-section of the column of Ainv
        for(k = 1; k < (Lp[i+1] - Lp[i]); k++){
          g = Lp[i] + k;
          // check if current index of L=0 (and hence not stored in CSC format)
          //// index beyond current row hence L=0
          //// first, catch up row within L to the non-zero indexes of L
          while( ((i + k) < Li[g]) && (g < Lp[i+1]) ) g++;
//Rprintf("\t\tkth row:%i | L row:%i | Ainv row:%i\n", i+k, Li[g], Ainvi[ Ainvp[j] + i + k ]);
          if( Li[g] == (i + k) ){ 
//FIXME: TEMPORARY checks/warnings
////  index beyond current column? (SHOULDN'T happen, given how for(k) setup)
if(g > Lp[i+1]){
  Rprintf("**** WARNING: `cs_chol2inv()` indexed beyond current column ****\n\n\n\n");
  continue;
} 

          raij += Lx[g] * Ainvx[ Ainvp[j] + i + k ];
        }
      }  // end for k

      // Ainv[i, j] = Divide `Sij - raij` by diagonal of L[i,i]
      Ainvij = (Sij - raij) / Lx[Lp[i]];
      // Add Ainvij to Ainv[i,j] and symmetric element Ainv[j, i]
      Ainvx[Ainvp[j]+i] = Ainvx[Ainvp[i]+j] = Ainvij;  
    }  // end for i
  }  // end for j

  Ainv->p = Ainvp; Ainv->i = Ainvi; Ainv->x = Ainvx;

 return(Ainv);

}
//TODO need a drop zeroes step?


*/
