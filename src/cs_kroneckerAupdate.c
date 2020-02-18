#include "gremlin.h"

void cs_kroneckerAupdate(const cs *G, const cs *A, const cs *C){

    int i, j, k, l, cnt,*Ap, an, gn;
    double  *Ax, *Gx;
	
    an = A->n ; Ap = A->p ;  Ax = A->x ;
    gn = G->n ; Gx = G->x ;

    // MEW 2020 02 17: cn and cnz set but not used
    // int cn, cnz, anz, gnz, 
    // cn = an*gn; 
    // anz = A->nzmax; gnz = G->nzmax; cnz = anz*gnz;
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



