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

  cs      *H, *Hinv;
  css     *sLh;
  csn     *Lh;

  csi     g, k, conP, si, si2, rw, boundFlag;

  double  d;

  if(!CS_CSC(A) || !wchBd) return (0);    // check arguments

  boundFlag = 0;
  conP = p;

  // Check for fixed (co)variance parameters
  for(k = 0; k < p; k++){
    if(con[k] == 0 || con[k] == 3){
      wchBd[k] = 1;  // only reset elements about to use
      conP--;
      if(con[k] == 3) boundFlag = 1;
    } else wchBd[k] = 0;
  }


  // make a generic working vector `w`
  // represents `grad` when boundFlag = 0 and `cdnu` when boundFlag = 1 
  double *w = new double[conP];

  if(boundFlag == 0){
    // remove fixed parameters from AI and dLdnu 
    //// (H for Hessian and w for gradient)
    si = 0;
    for(k = 0; k < p; k++){
      if(con[k] == 0) continue;
        w[si] = dLdnu[k];
        newnu[si] = nu[k];
        si++; 
    }
  }
  H = cs_droprowcol(A, wchBd);      





  // Check if Hessian can be inverted
  sLh = cs_schol(1, H);
  Lh = cs_chol(H, sLh);
  if(Lh == NULL){
    f = 3e-5; 

    if(v > 1){
      Rprintf("\n\tH cholesky decomposition failed:\n\t   Hessian matrix may be singular - modifying diagonals and re-trying");
        if(v > 2) Rprintf("\n\tH modification: %6.3g\n", f);
    } // end if v>1
        
  }  // end if Lh=NULL



  // Check/modify H matrix to 'ensure' positive definiteness
  /* `f` is factor to adjust H matrix
     (e.g., based off Meyer 1997 eqn 58 and WOMBAT manual A.5 strategy 3b)
     HOWEVER, differs from above in that I do not use eigenvalues here,
     rather the diagonals of the Cholesky factorization are inspected
       see `remlIt.gremlinR()` for approach like Meyer (with eigenvalues) */
/*
See `La_rs` defined around line 153 of "R/src/modules/lapack/Lapack.c"
What R's `eigen()` calls
 */

  //// check if any small/zero diagonals on Cholesky of H (Lh)
  ////// If any diagonals of Lh < 0 modify H by adding f to diagonals
  if(Lh != NULL){
    for(g = 0; g < Lh->L->n; g++){
      d = Lh->L->x[ Lh->L->p[g] ];
      if(d < ezero[0]){
        f = 3e-5;
        if(v > 2) Rprintf("\n\tSmall diagonal on Cholesky of H: adding %6.3g\n", f);

      }
    }
  }
  ////XXX ASSUME H is full matrix so H->p[g] + g = diagonal
  for(g = 0; g < H->n; g++) H->x[ H->p[g] + g ] += f;

  // Now check AGAIN if Hessian can be inverted
  cs_sfree(sLh); 
  cs_nfree(Lh);
  sLh = cs_schol(1, H);
  Lh = cs_chol(H, sLh);
  if(Lh == NULL){
    if(v > 1){
      Rprintf("\n\tH cholesky decomposition failed:\n\t   Hessian matrix may be singular");
    } // end if v>1

    cs_sfree(sLh);
    cs_nfree(Lh);
    cs_spfree(H);
    delete [] w;
   return (0);

  } //<-- end if H cannot be inverted  


  Hinv = cs_inv(H);
  if(Hinv == NULL){
    cs_sfree(sLh); 
    cs_nfree(Lh);
    cs_spfree(H);
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
      for(k = 0; k < H->m; k++) w[k] = 0.0;
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
      cs_sfree(sLh); 
      cs_nfree(Lh);
      cs_spfree(Hinv);
      cs_spfree(H);

 return (1) ;	/* success; return 1 */
}



