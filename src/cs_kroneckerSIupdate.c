#include "MCMCglmm.h"

void cs_kroneckerSIupdate(const cs *A, int nI, const cs *C){

    int i, j, k, cnt;

    cnt = 0;	
    for (i = 0 ; i <A->n; i++){
      for (j = 0 ; j < nI ; j++){
        for(k=A->p[i];  k < A->p[i+1] ; k++){
          C->x[cnt] = A->x[k];
          cnt++;
        } 
      }     
    }

}



