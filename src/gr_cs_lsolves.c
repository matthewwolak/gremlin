#include "gremlin.h"
/*********************************************************/
/*  MEW ADAPTED from below: suit special case in gremlin */
/** Automatically skips over zeroes in `x`              **/
/*  CSparse source code: Sept 12, 2017: version 3.2.0	 */
/** from SuiteSparse-5.1.0				**/
/** See License in ./src/CSparse/Doc/License.txt	**/
/*********************************************************/


/* solve Ax=k where Lx=b, L'b=k, and x and b are dense.
   x=b on input, solution on output. */
// combines cs_lsolve followed by cs_ltsolve
//// start at element k of x
csi gr_cs_lltsolve (const cs *L, double *x, csi k)
{
    csi p, j, n, *Lp, *Li ;
    double *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;

    // lsolve
    for (j = k ; j < n ; j++){
        if(x[j] != 0.0){
          x [j] /= Lx [Lp [j]] ;
          for (p = Lp [j]+1 ; p < Lp [j+1] ; p++) {
            x [Li [p]] -= Lx [p] * x [j] ;
          }
        }  // end if NOT 0
    }


    // ltsolve
    for (j = n-1 ; j >= 0 ; j--){
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [j] -= Lx [p] * x [Li [p]] ;
        }
        x [j] /= Lx [Lp [j]] ;
    }
    return (1) ;
}


