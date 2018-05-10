#include "MCMCglmm.h"

cs *cs_kroneckerI(const cs *A, int nI){

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
	
    for (j = 0 ; j < cn ; j++){
         for (i = 0 ; i < am ; i++){
            Ci[cnt] = i*nI+j%nI;
            cnt++;
         }      
    }

    cnt = 0;
    Cp[0] = 0;
    for(i = 0; i<an; i++){
     	for (j = 0 ; j < nI ; j++){
          cnt++;
          Cp[cnt] = Cp[cnt-1]+am;
	}
    }

    cnt = 0;
    for(i = 0; i < an; i++){
     	for(j = 0 ; j < nI ; j++){
            for(k = 0; k < am; k++){
              Cx[cnt] = Ax[i*an+k];
              cnt++;
            }
	}
    }

    cs_sprealloc (C, 0) ;		// remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



