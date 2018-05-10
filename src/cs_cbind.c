#include "MCMCglmm.h"

cs *cs_cbind(const cs *A, const cs *B){

    int j, anz, bnz, cnz, *Cp, *Ap, *Bp, *Ci, *Ai, *Bi, cm, an, cn;
    double *Ax, *Bx, *Cx;
    cs *C;
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;	 /* check inputs */
    if (A->m != B->m) return (NULL) ;
    cm = A->m ; an = A->n; cn = (A->n)+(B->n);
    anz = A-> nzmax; bnz = B-> nzmax; cnz = anz+bnz;
    Ap = A->p; Bp = B->p; 
    Ai = A->i; Bi = B->i; 
    Ax = A->x; Bx = B->x; 
    C = cs_spalloc (cm, cn, cnz, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0)) ;
//    if (!C ) error("cs_cbind out of memory"); 
    Ci = C->i;
    Cx = C->x;
    Cp = C->p;
    for (j = 0 ; j < an ; j++){Cp [j] = Ap[j];}
    for (j = an ; j < cn ; j++){Cp [j] = Bp[j-an]+anz;}
    for (j = 0 ; j < anz ; j++){
    	Ci [j] = Ai[j];
   	Cx [j] = Ax[j];
    }
    for (j = anz ; j < cnz ; j++){
    	Ci [j] = Bi[j-anz];
   	Cx [j] = Bx[j-anz];
    }
    Cp [cn] = cnz ;			//finalize the last column of C
    cs_sprealloc (C, 0) ;		// remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



