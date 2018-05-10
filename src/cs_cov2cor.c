#include "MCMCglmm.h"

void cs_cov2cor(const cs *A){
    
    int m, j,k;
    m = A->n;
	for(j=0; j<m; j++){
		for(k=0; k<m; k++){
			if(j!=k){
			  A->x[j*m+k] = A->x[j*m+k]/sqrt(A->x[j*m+j]*A->x[k*m+k]);
		    }
		}
	}
    for(j=0; j<m; j++){
		A->x[j*m+j] = 1.0;
	}		
}


