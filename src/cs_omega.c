#include "MCMCglmm.h"

cs *cs_omega(cs **KGinv, int nG, cs *pvB){

    int i, j, *Cp, *Ap, *Ci, *Ai, an, anz, cn=0, cnz=0, nlG=0, nlE=0, nlB=0, cnt=0, cnt2=0;
    double *Cx, *Ax;
    cs *C;

    nlB = pvB->n;  // Note nlB is the length of the fixed effect vetcor NOT fixed effects within traits

    for(i = 0; i<nG; i++){
      cn  += KGinv[i]->n;
      cnz += KGinv[i]->nzmax;
    }

    cn += nlB;
    cnz+= pvB->nzmax;

    C = cs_spalloc (cn, cn, cnz, 1, 0) ;	 /* allocate result */

    if (!C ) return (cs_done (C, NULL, NULL, 0));  
//      if (!C) error("cs_omega out of memory"); 

    Cp = C->p ; Ci = C->i ; Cx = C->x;
   
    for (j = 0 ; j<nlB; j++){ 
       Cp[cnt] = cnt2+pvB->p[j];
       cnt++;
    }
    cnt2+=pvB->nzmax;

    for (j = 0 ; j < pvB->nzmax; j++){
       Ci[j] = pvB->i[j];
       Cx[j] = pvB->x[j];
    }      
    nlG += nlB;
    nlE += pvB->nzmax;

     for(i = 0; i<nG; i++){

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



