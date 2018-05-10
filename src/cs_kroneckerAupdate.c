#include "MCMCglmm.h"

void cs_kroneckerAupdate(const cs *G, const cs *A, const cs *C){

    int i, j, k, l, cnt, anz, gnz, cnz, *Ap, an, gn, cn;
    double  *Ax, *Gx;
	
    an = A->n ; anz = A->nzmax; Ap = A->p ;  Ax = A->x ;
    gn = G->n ; gnz = G->nzmax; Gx = G->x ;

    cn = an*gn; cnz = anz*gnz;
    cnt = 0;	
    for (i = 0 ; i < gn ; i++){
      for (j = 0 ; j < an ; j++){    
        for (k = 0 ; k < gn ; k++){    
          for (l = Ap[j] ; l < Ap[j+1] ; l++){
            C->x[cnt] = Ax[l]*Gx[i*gn+k];
            cnt++;
          }
        }
      }
    }      
}



