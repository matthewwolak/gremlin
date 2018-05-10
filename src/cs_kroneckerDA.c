#include "MCMCglmm.h"

cs *cs_kroneckerDA(double *D, int n, const cs *A){

    int i, j, cnt, anz, an;

    cs *C;
    if (!CS_CSC (A)) return (NULL);                         
    an = A->n ; anz = A->p[an]; 
    C = cs_spalloc (an*n, an*n, anz*n, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));  

    cnt = 0;	
    for (i = 0 ; i < n ; i++){
      for (j = 0 ; j < anz ; j++){    
        C->x[cnt] = A->x[j]*D[i];
        C->i[cnt] = A->i[j]+i*an; 
        cnt++;
      }
    }
    cnt=0;
    for (i = 0 ; i < n ; i++){
      for (j = 0 ; j < an; j++){ 
        C->p[cnt] = A->p[j]+i*(A->p[an]); 
        cnt++;
      }
    }      
    C->p[n*an] = n*anz;

    cs_sprealloc (C, 0);		        // remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



