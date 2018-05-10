#include "MCMCglmm.h"

cs *cs_kroneckerA(const cs *G, const cs *A){

    int i, j, k, l, cnt, anz, gnz, cnz, *Cp, *Ap, *Gp, *Ci, *Ai, *Gi, an, gn, cn;
    double *Cx, *Ax, *Gx;
    cs *C;
    if (!CS_CSC (A)) return (NULL);                         
    an = A->n ; anz = A->nzmax; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    gn = G->n ; gnz = G->nzmax; Gp = G->p ; Gi = G->i ; Gx = G->x ;

    cn = an*gn; cnz = anz*gnz;
    C = cs_spalloc (cn, cn, cnz, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));  
    Cp = C->p ; Ci = C->i ; Cx = C->x ;   
    cnt = 0;	
    for (i = 0 ; i < gn ; i++){
      for (j = 0 ; j < an ; j++){    
        Cp[j+an*i] = gn*(Ap[j]+anz*i);
        for (k = 0 ; k < gn ; k++){    
          for (l = Ap[j] ; l < Ap[j+1] ; l++){
            Ci[cnt] = Ai[l]+k*an;
            Cx[cnt] = Ax[l]*Gx[i*gn+k];
            cnt++;
          }
        }
      }
    }      
    Cp[cn] = cnz;
    cs_sprealloc (C, 0) ;		        // remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



