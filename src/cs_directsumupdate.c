#include "MCMCglmm.h"

void cs_directsumupdate(cs **KGinv, int nG, int nGR, const cs *C){

    int i, j, an, anz, nlG=0, nlE=0;
    double *Ax;

	for(i = nG; i<nGR; i++){

        an = KGinv[i]->n ; anz = KGinv[i]->nzmax; Ax = KGinv[i]->x;

        for (j = 0 ; j < anz; j++){
          C->x[j+nlE] = Ax[j];
        }      
        nlG += an;
        nlE += anz;
    }
}



