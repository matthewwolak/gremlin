#include "MCMCglmm.h"

void cs_kroneckerIupdate(const cs *A, int nI, const cs *C){

    int i, j, k, cnt, an, am;
    double *Ax;

    an = A->n;
    am = A->m;
    Ax = A->x;
	
    cnt = 0;
    for(i = 0; i < an; i++){
     	for(j = 0 ; j < nI ; j++){
            for(k = 0; k < am; k++){
              C->x[cnt] = Ax[i*an+k];
              cnt++;
            }
	}
    }
}



