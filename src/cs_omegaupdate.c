#include "MCMCglmm.h"

void cs_omegaupdate(cs **KGinv, int nG, cs *pvB, const cs *C){

    int i, j, anz, nlE=0;
    double  *Ax;
     
    nlE += pvB->nzmax;

     for(i = 0; i<nG; i++){
		 
         anz = KGinv[i]->nzmax; Ax = KGinv[i]->x;

        for (j = 0 ; j < anz; j++){
          C->x[j+nlE] = Ax[j];
        }      
        nlE += anz;
    }
}



