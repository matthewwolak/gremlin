#include "MCMCglmm.h"

cs *cs_kroneckerDI(double *D, int n, int nI){

    int i, j, cnt;
    cs *C;

    C = cs_spalloc (n*nI, n*nI, n*nI, 1, 0) ;	 /* allocate result */
    if (!C ) return (cs_done (C, NULL, NULL, 0));   

    cnt=0;
    for(i = 0; i<n; i++){	
      for (j = 0 ; j < nI ; j++){
         C->i[cnt] = cnt;
         C->p[cnt] = cnt;
         C->x[cnt] = D[i];
         cnt++;
      }
    }
    C->p[n*nI] = n*nI;
    cs_sprealloc (C, 0) ;		        // remove extra space from C 
    return (cs_done (C, NULL, NULL, 1)) ;	/* success; free workspace, return C */
}



