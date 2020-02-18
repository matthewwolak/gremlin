#include "gremlin.h"




/* AI returned */
cs *cs_ai(const cs *BLUXs, double *nu,
        const cs *R, const cs *KRinv, const cs *tWKRinv,
        double *rory,  // residuals if lambda=FALSE else y if lambda=TRUE
        const cs *W, const cs *tW, csi n, csi p, csi nG, csi *rfxlvls, csi nb,
	const cs *Lc, const csi *Pinv,
	csi thetaR,       // 0 if lambda=TRUE
        double sigma2e,    // 1.0 if lambda=FALSE
	double ezero
){

  int     lambda;
  double  sln_k;
  cs      *AI, *Rinv, *B, *tB, *BRHS, *tBRinvB, *tBKRinv, *tS, *Scol, *tSBRHS;
  csi     g, i, k, cnt, si, qi, ei;

  if(!CS_CSC(BLUXs) || !nu) return(0);

  double  *tmp_sln = new double[BLUXs->m];
  if(!tmp_sln) return(0);

  if(thetaR != 0 && fabs(sigma2e - 1.00) < ezero) lambda = 0; else lambda = 1;

  if(lambda == 1){
    Rinv = cs_spalloc(1, 1, 1, true, false);
    Rinv->i[0] = 0; Rinv->p[0] = 0; Rinv->p[1] = 1; Rinv->x[0] = (1.0 / sigma2e);
  }else{
    Rinv = cs_inv(R);
  }
  B = cs_spalloc(n, p, (n * p), true, false);
    B->p[0] = 0;

  si = nb;

  for(g = 0; g < nG; g++){
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    // Meyer 1997: eqn. 20 [also f(theta) in Johnson & Thompson 1995 eqn 9c)]
    // Knight 2008: eqn 2.32-2.35 and 3.11-3.13 (p.34)
    //TODO for covariances see Johnson & Thompson 1995 eqn 11b

    //// Go column-by-column of part of W corresponding to gth Z
    ////// assign inverse-(co)variance * kth solution to each individual
//FIXME TODO: think currently ASSUMES each Z has only ONE non-zero per row
//FIXME
//FIXME just send list of Zs to c++ then this function and do actual matrix multiplications
    ////// B[,g] = Z_g %*% sln[si:ei, 1] %*% Ginv[g] FIXME is order correct? See difference in order between Johnson & Thompson (e.g., eqn. 11b) and Meyer 1997 eqn 20
    cnt = 0;
    for(k = si; k < ei; k++){
      sln_k = BLUXs->x[k];
      for(i = W->p[k]; i < W->p[k+1]; i++){
        B->i[cnt] = W->i[i];
        B->x[cnt] = W->x[i] * sln_k * (1.0 / nu[g]);
        cnt++;
      }  // end for i (each row in a column of Z - share same solution) 
      B->p[k-si+1] = cnt;
    }  // end for k (each solution for a given gth parameter)
    si = ei + 1;
  }  // end for g

  //FIXME TODO Check what to do if more than 1 residual variance parameter
  for(k = 0; k < n; k++){
    B->i[cnt] = k;
    B->x[cnt] = rory[k] * Rinv->x[0];  //TODO ?? what to do if Rinv is a matrix?
    cnt++;
  }
  B->p[p] = cnt;
  tB = cs_transpose(B, true);
  
  // Set up modified MME like the MMA of Meyer 1997 eqn. 23
  //// Substitute `B` instead of `y`

  // lambda: (co)variance ratios
  if(lambda == 1){
    BRHS = cs_multiply(tW, B);
    tBRinvB = cs_multiply(tB, B);

  //// NOT lambda: Straight (co)variances
  }else{
    tBKRinv = cs_multiply(tB, KRinv);
    // next is actually `tBKRinvB`, want 1 name for this and when lambda=TRUE
    tBRinvB = cs_multiply(tBKRinv, B);
      cs_spfree(tBKRinv);
    BRHS = cs_multiply(tWKRinv, B);
  }    // end if lambda=FALSE

  // tBPB
  //// Johnson & Thompson 1995 eqns 8,9b,9c
  ////// (accomplishes same thing as Meyer 1997 eqns 22-23 for a Cholesky
  ////// factorization as Boldman & Van Vleck eqn 6 applied to AI calculation)
  // directly create transpose to pre-multiply with BRHS in `... - crossprod()`
  tS = cs_spalloc(BRHS->n, BRHS->m, BRHS->p[BRHS->n], true, false);

  // create temporary "row" of tS as a 1 column matrix
  Scol = cs_spalloc(BRHS->m, 1, BRHS->m, true, false);
  Scol->p[0] = 0;

  cnt = 0;
  for(k = 0; k < BRHS->m; k++){
    tS->p[k] = cnt;
    for(i = 0; i < p; i++){
      tS->i[cnt] = i;
      tS->x[cnt] = 0.0;
      cnt++;
    }  // end rows in column k
    // Now deal with Scol (just 1 column)
    Scol->i[k] = k;
    Scol->x[k] = 0.0;
  }
  tS->p[tS->n] = cnt;
  Scol->p[1] = Scol->m;

  for(k = 0; k < p; k++){
    cnt = k * BRHS->m;
    for(i = 0; i < BRHS->m; i++){
      tmp_sln[k] = BRHS->x[cnt+i];
    }

    cs_ipvec(Pinv, tmp_sln, Scol->x, BRHS->m);	     // x = P*b 
    cs_lsolve(Lc, Scol->x);                          // x = L\x 
    cs_ltsolve(Lc, Scol->x);		             // x = L'\x 
    cs_pvec(Pinv, Scol->x, tmp_sln, BRHS->m);        // b = P'*x 
  
    // put tmp_sln into kth row of tS
    for(i = 0; i < BRHS->m; i++){
      tS->x[i*p + k] = tmp_sln[i];
    }  
  }

  tSBRHS = cs_multiply(tS, BRHS);
  // AI = -0.5 * (tBRinvB - tS %*% BRHS)
  //// `alpha*A + beta*B` gives 1*tBRinvB + -1*tSBRHS
  AI = cs_add(tBRinvB, tSBRHS, 1.0, -1.0);  
  for(i = 0; i < AI->p[AI->n]; i++) AI->x[i] *= 0.5;
  if(lambda == 1) for(i = 0; i < AI->p[AI->n]; i++) AI->x[i] /= sigma2e;


  cs_spfree(Rinv);
  cs_spfree(B);
  cs_spfree(tB);
  cs_spfree(BRHS);
  cs_spfree(tBRinvB);  
  cs_spfree(tS);
  cs_spfree(Scol);
  cs_spfree(tSBRHS);

//  delete [] tmp_sln;
// return(AI);
  // success, free tmp_sln, return AI and 1=success
 return(cs_done(AI, tmp_sln, NULL, 1));
}



