#include "MCMCglmm.h"

cs *cs_kroneckerSI(const cs *A, int nI){

/***************************************************************************************************/
/* Matrices are stored in compressed-column format where i indicates the row indicies of non-zero  */
/* values, p indexes the first i of each column, and x the actual non-zero values.  For example,if */
/*                                                                                                 */
/*       0 1 0            p = {0, 1, 3, 4}                                                         */
/*   X = 1 0 2    then,   i = {1, 0, 2, 1}                                                         */
/*       0 1 0            x = {1.0, 1.0, 1.0, 2.0}                                                 */
/*                                                                                                 */
/* where the final element of p is always length(i).                                               */
/* dim is a vector with the number of rows and columns, and nzmax the max number of non-zero values*/
/***************************************************************************************************/


    int i, j, k, cnt, anz, cnz, *Cp, *Ap, *Ci, *Ai, am, an, cm, cn;
    double *Cx, *Ax;
    cs *C;
    if (!CS_CSC (A)) return (NULL);                         
    am = A->m ; an = A->n ; anz = A->nzmax; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    cm = am*nI; cn = an*nI; cnz = anz*nI;
    C = cs_spalloc (cm, cn, cnz, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));   

    Cp = C->p ; Ci = C->i ; Cx = C->x ;   
    cnt = 0;	
    for (i = 0 ; i <an; i++){
      for (j = 0 ; j < nI ; j++){
        for(k=A->p[i];  k < A->p[i+1] ; k++){
          Ci[cnt] = A->i[k]*nI+j;
          Cx[cnt] = A->x[k];
          cnt++;
        } 
      }     
    }
    cnt = 0;
    Cp[0] = 0;
    for(i = 1; i<=an; i++){
     	for (j = 0 ; j < nI ; j++){
          cnt++;
          Cp[cnt] = Cp[cnt-1]-A->p[i-1]+A->p[i];
	}
    }

    cs_sprealloc (C, 0) ;		// remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */}



