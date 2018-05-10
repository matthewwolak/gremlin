#include "MCMCglmm.h"

cs *cs_directsum(cs **KGinv, int nG, int nGR){

    int i, j, *Cp, *Ap, *Ci, *Ai, an, anz, cn=0, cnz=0, nlG=0, nlE=0, cnt=0, cnt2=0;
    double *Cx, *Ax;
    cs *C;

    for(i = nG; i<nGR; i++){
      cn  += KGinv[i]->n;
      cnz += KGinv[i]->nzmax;
    }

    C = cs_spalloc (cn, cn, cnz, 1, 0) ;	 /* allocate result */

    if (!C ) return (cs_done (C, NULL, NULL, 0));   
//     if (!C ) error("cs_directsum out of memory");  

    Cp = C->p ; Ci = C->i ; Cx = C->x;
   
   for(i = nG; i<nGR; i++){

        if (!CS_CSC (KGinv[i])) return (NULL); 
                        
        an = KGinv[i]->n ; anz = KGinv[i]->nzmax; Ap = KGinv[i]->p; Ai = KGinv[i]->i; Ax = KGinv[i]->x;

        for (j = 0 ; j < an; j++){
          Cp[cnt] = cnt2+Ap[j];
          cnt++;
        } 
        cnt2+=anz;

        for (j = 0 ; j < anz; j++){
          Ci[j+nlE] = Ai[j]+nlG;
          Cx[j+nlE] = Ax[j];
        }      
        nlG += an;
        nlE += anz;
    }
    Cp[cn]=cnz;
    cs_sprealloc (C, 0) ;		// remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



