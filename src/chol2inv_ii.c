#include "gremlin.h"


/* sparse matrix inversion from Cholesky Decomposition: L t(L)=A  */
/* Follows the method/equations by:
  Takahashi, Fagan, & Chin. 1973. Formation of a sparse bus impedance matrix
  and its application to short circuit study. 8th PICA Conference Proceedings,
  Minneapolis, MN.

and implementation in SuiteSparse --> MATLAB_Tools/sparseinv by Tim Davis.
SuiteSparse:::sparseinv works on L D L' = A Cholesky factorization, whereas, below the code is adapted to work from L L' factorization.
 
Target is to only keep/return diagonal elements of inverse if L L' = A, this
  gives (diag(Ainv) by only calculating elements of Ainv based on non-zero
  elements of L/L'
*/

/* replaces Cinv_ii with diagonals from partially inverted matrix  */
/* requires non-zero pattern of L + L' already known in matrix Z   */

csi chol2inv_ii(const cs *L, const csi *Pinv,
	cs *Z, int *Zdiagp, double *Cinv_ii, int itErr){  

  csi     i, j, k, p, up, zp, Ln, *Lp, *Li, *Zp, *Zi, *Lmunch,
          znz;
  double  zkj, ljk, *Lx, *Zx, *d, *z;

  Ln = L->n;
  Lp = L->p; Li = L->i; Lx = L->x;
  Zp = Z->p; Zi = Z->i; Zx = Z->x;
  
  d = cs_malloc(Ln, sizeof (double));     // diagonal of L
  z = cs_malloc(Ln, sizeof (double));    // holds Z[, j] workspace
  Lmunch = cs_calloc(Ln, sizeof (csi));  // point to last (working) entry L[, k]


  // If not already initialized, fill Z/prtCinv with non-zero pattern L+L'
  if(itErr == 0){ 
    // first determine Zp (briefly what cs_compress does at beginning)
    //// start with column counts of L BUT subtract 1 for diagonal
    ////// use Lmunch to keep column counts
    for(k = 1; k <= Ln; k++) Lmunch[k-1] = Lp[k] - Lp[k-1] - 1;
    // to figure column counts of L' add kth row from Li to kth column count
    //// will add to count for diagonal (why subtracted 1 above)
    for(k = 0; k < Lp[Ln]; k++) Lmunch[ Li[k] ]++ ;  // end with column counts
    cs_cumsum(Zp, Lmunch, Ln);
    znz = Zp[Ln];  // 1+total number non-zeroes in Z


    // fill upper triangle of Z with L' (though don't make L')
    //// L[j, k] is L'[k, j]/Z[k, j] so add this to jth column of Z
    //// starting from top row means will place Zi in correct row order
    //// use Lmunch to keep track of next Zp[k] pointer to use for each column
    for(k = 0; k < Ln; k++) Lmunch[k] = Zp[k];
    for(k = 0; k < Ln; k++){
      Zdiagp[k] = Lmunch[k];// current Lmunch should=Zp location for Z[k, k]
      for(j = Lp[k]; j < Lp[k + 1]; j++){
        i = Li[j];         // L[i, k] means this goes to Z[k, i]
        Zi[ Lmunch[i]++ ] = k;
      }  // end for j
    }  // end for k
        

    // fill lower triangle of Z with L; column-by-column
    for(k = 0; k < Ln; k++){
      // skip diagonal (start at Lp[k] + 1)
      for(j = Lp[k] + 1; j < Lp[k + 1]; j++) Zi[ Lmunch[k]++ ] = Li[j];
    }  // end for k
      
    // check that non-zero diagonals exist
    for(k = 0; k < Ln; k++){
      // will ID Z_jj/diagonal as zero=not allowed
      if(Zi[ Zdiagp[k] ] != k) itErr = -1;
    }  // end for k

  }  // end itErr=0 and need to initialize Zi and Zp




  if(itErr >= 0){
    znz = Zp[Ln];  // 1+total number non-zeroes in Z (even if already done above)
    for(p = 0; p < znz; p++) Zx[p] = 0.0;  // set Zx to zero/clear previous values
  
    // initialize Z diagonals and d & fill Lmunch[j] with last entry in L[, j]
    for(j = 0; j < Ln; j++){
      zkj = 1.0 / Lx[ Lp[j] ];
      d[j] = zkj;
      Zx[ Zdiagp[j] ] = zkj * zkj;  // using LL' (if LDL' then zkj/d[j]) 
      Lmunch[j] = Lp[j + 1] - 1;
    }  // end for j



    // Compute sparse inverse subset
    for(j = Ln-1; j >= 0; j--){
  
      // scatter Z[, j] into z workspace (TODO worth using cs_scatter?)
      //// only lower triangular part is needed, since upper is all zero to start
      for(p = Zdiagp[j]; p < Zp[j+1]; p++) z[ Zi[p] ] = Zx[p];


      // compute Z[, j]  ***above*** diagonals
      //// for ROW k=(j-1) to 0 but only for non-zero entries Z[k, j]
      //// moves from just above diagonal up to top of row of Z[, j]
      for(p = Zdiagp[j]-1; p >= Zp[j]; p--){

        // Z[k, j] = -U[k, (k+1):n] * Z[(k+1):n, j]  
        k = Zi[p];
        zkj = 0.0;
      
        /* U=L' and if U stored by row then same as L stored by column (so L=U)
          at Z[k, j] starts in k row of U (col of L) & moves to end of U(L)   */
        for(up = Lp[k]; up < Lp[k + 1]; up++){
          i = Li[up];  // skip diagonal of U 
          if(i > k){
            // need to adjust L (LL') so it is L from LDL' by "... * (1/D[k])"
            zkj -= Lx[up] * d[k] * z[i];
          }
      
        }  // end for up
        z[k] = zkj;
      
      }  // end for p
    
    
    
      // Left-looking update to Z ***below*** diagonals
      //// for ROW k=(j-1) to 0 but only for non-zero entries Z[k, j]
      //// moves from just above diagonal up to top row of Z[, j]
      for(p = Zdiagp[j]-1; p >= Zp[j]; p--){
    
        k = Zi[p];      // upper triangle row of Z
        // ljk = L[j, k]
        if((Lmunch[k] < Lp[k]) || (Li[ Lmunch[k] ] != j)){
          continue;  // skip L[j, k] = 0 since no math/work to do
        }
      
        // need to adjust L (LL') so it is L from LDL' by "... * (1/D[k])"
        ljk = d[k] * Lx[ Lmunch[k]-- ];
      
        // Zp(k+1):n, k] = Z[(k+1):n, k] - Z[(k+1):n, j] * L[j, k]
        //// moves from diagonal of Z[, k] to last row of Z[, k]
        for(zp = Zdiagp[k]; zp < Zp[k+1]; zp++){ 
          Zx[zp] -= z[ Zi[zp] ] * ljk;
        }  // end for zp
      
      }  // end for p



      // Gather Z[, j] back from z workspace 
      //// moves down column Z[, j] by (non-zero) rows
      for(p = Zp[j]; p < Zp[j+1]; p++){
        i = Zi[p];
        Zx[p] = z[i];
        z[i] = 0.0;
      }  // end for p
      
    }  // end for j
      

  
    // Fill diagonals of Z into vector Cinv_ii
    //// Z is in order of sLc/Lc, but Cinv should be original ordering
    for(j = 0; j < Ln; j++){
      k = Zdiagp[ Pinv[j] ];       /* essentially `cs_pvec` */
      Cinv_ii[j] = Zx[k];
    }
  
  }  // end itErr >= 0 (no errors/zero diagonals)
  
  cs_free(d);
  cs_free(z);
  cs_free(Lmunch);

 if(itErr >= 0){
   return(1);
 } else return(0); 

}






