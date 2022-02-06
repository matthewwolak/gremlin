#include "gremlin.h"

// Find next set of parameters using a quasi-Newton method/algorithm
//// Meyer 1989 pp. 326-327 describes quasi-Newton methods 
// see Meyer 1997 eqn 58 for Marquardt 1963: theta_t+1=theta_t - (H_t + k_t * I)^{-1} g_t 
// What I do below is similar: except k_t=f
//// Mrode 2005 eqn 11.4
//// Johnson and Thompson 1995 eqn 12
//////(though gremlin uses `+` instead of J & T '95 `-` because
////// gremlin multiplies gradient by -0.5 in `gradFun()`)


/* 1 returned if successful, else NULL
   `newnu` replaced with next set of nu parameters   */
csi qNewtRhap(double *nu, double *newnu, double *dLdnu, const cs *A,
	      csi p, csi *con, csi *wchBd, double f, double *ezero, csi v)
{

  cs      *Hinv;
  csi     g, k, conP, si, si2, rw, boundFlag;

  if(!CS_CSC(A) || !wchBd) return (0);    // check arguments

  boundFlag = 0;
  conP = p;

  // Check for fixed (co)variance parameters
  for(k = 0; k < p; k++){
    if(con[k] == 0 || con[k] == 3){
      conP--;
      if(con[k] == 3) boundFlag = 1;
    }
  }


  // make a generic working vector `w`
  // represents `grad` when boundFlag = 0 and `cdnu` when boundFlag = 1 
  double *w = new double[conP];

  if(boundFlag == 0){
    // remove fixed parameters from AI and dLdnu 
    //// (w for gradient)
    si = 0;
    for(k = 0; k < p; k++){
      if(con[k] == 0) continue;
        w[si] = dLdnu[k];
        newnu[si] = nu[k];
        si++; 
    }
  }


  Hinv = cs_inv_withDiagMod(A, con, wchBd, ezero, v);
  if(Hinv == NULL){
    delete [] w;
    return (0);
  }


  // fill `newnu` with parameters proposed for next iteration
  if(boundFlag == 0){
    ////  cs_gaxpy is y = A*x+y where y=newnu and x=w/gradient
    cs_gaxpy(Hinv, w, newnu);



  } else{

      /* matrix multiplications "by hand" to pull out subsets
         newnu[-bad] += invH_uu %*% (grad_u - H_uc %*% (newnu[bad]-nu[bad]))  */
                 
      // `w` = conditional `dnu`
      for(k = 0; k < Hinv->m; k++) w[k] = 0.0;
      // H_uc %*% (newnu[bad] - nu[bad])
      //// in practice: w += H_uc[, k] * (newnu[bad] - nu[bad])[k]
      for(g = 0; g < p; g++){  // go through original AI
        if(con[g] != 3) continue;  // only g for a boundary parameter
        si = g;  // used to index gth parameter of AI in newnu and w
        for(rw = 0; rw < g; rw++) if(con[rw] == 0) si--;
          si2 = 0;  // used to index row of H_uc[, k] in w
          for(k = A->p[g]; k < A->p[g+1]; k++){
            if(wchBd[ A->i[k] ] == 0){
              w[si2] += A->x[k] * (newnu[si] - nu[g]);
              si2++;
            }
          }  //<-- end for k (rows of AI column (g) for boundary parameter) 
        }  //<-- end for g (columns of AI) 
            
        // ... (grad_u - ...) operation
        // at same time replace newnu[-c(fixed,bad)] with nu[-c(fixed,bad)]
        //// w has elements 0-to-(No. non-fixed & non-boundary)
        //// newnu has elements 0-to-(No. non-fixed)
        si = 0; si2 = 0;
        for(g = 0; g < p; g++){
          if(con[g] == 0) continue;  // assumes either fixed or boundary
          if(con[g] != 3){
            w[si] = dLdnu[g] - w[si];
            newnu[si2] = nu[g];  // don't replace boundary (already set above)
            si++;
          }
          si2++;
        }

        // invH_uu %*% (...)
        rw = 0;  // for indexing the current row of invH_uu to work across
        si = 0;  // index newnu (only contiains non-fixed parameters)
        for(g = 0; g < p; g++){
          if(con[g] == 0) continue;
          if(con[g] != 3){
          // go across columns for row `rw`
            for(k = 0; k < Hinv->n; k++){
              if(Hinv->i[ Hinv->p[k] + rw ] == rw){
                newnu[si] += Hinv->x[ Hinv->p[k] + rw] * w[k];
              }
            }  // end for k
            rw++;
          }  // end if NOT a boundary
          si++;
        }  // end for g 

      }  // end if/else boundaryFlag

     

      // Cleanup:
      delete [] w;
      cs_spfree(Hinv);

 return (1) ;	/* success; return 1 */
}



