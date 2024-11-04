#include "gremlin.h"


/* replaces AI */
cs *ai(cs *BLUXs, cs **Ginv,
        cs *R, cs *KRinv, cs *tWKRinv,
        double *rory,  // residuals if lambda=FALSE else y if lambda=TRUE
        cs *W, cs *tW, csi n, csi p, csi nG, csi *rfxlvls, csi nb,
	cs *Lc, csi *Pinv,
	csi thetaR,       // 0 if lambda=TRUE
        double sigma2e    // 1.0 if lambda=FALSE
){

  int     lambda;
  double  *BRHSx;
  cs      *AI, *Rinv, *Zg, *B, *tB, *BRHS, *tBRinvB, *tBKRinv,
	  *S, *tS, *tSBRHS;
  csi     g, i, k, cnt, si, qi, ei, *Sp;

  if(!CS_CSC(BLUXs)) return (0);

  double  *p_sln = new double[BLUXs->m];
  double  *Scol = new double[BLUXs->m];
  if(!p_sln || !Scol){
    delete [] p_sln;
    delete [] Scol;
    return (0);
  }

  if(thetaR != 0 && fabs(sigma2e - 1.00) < DBL_EPSILON) lambda = 0; else lambda = 1;

  if(lambda == 1){
    Rinv = cs_spalloc(1, 1, 1, true, false);
    Rinv->i[0] = 0; Rinv->p[0] = 0; Rinv->p[1] = 1; Rinv->x[0] = (1.0 / sigma2e);
  }else{
    Rinv = cs_inv(R);
    if(Rinv == NULL){
      delete [] p_sln;
      delete [] Scol;
      return (0);
    } 
  }
  B = cs_spalloc(n, p, (n * p), true, false);
  // Fill B with all zeroes
  cnt = 0;
  B->p[0] = 0;
  for(k = 0; k < p; k++){
    for(g = 0; g < n; g++){
      B->i[n*k + g] = g;     
      B->x[n*k + g] = 0.0;
      cnt++;
    }
    B->p[k+1] = cnt;
  }

  si = nb;

  for(g = 0; g < nG; g++){
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    Zg = cs_spalloc(n, qi, W->p[ei+1] - W->p[si], true, false);

    cnt = 0;
    for(k = W->p[si]; k < W->p[ei+1]; k++){
      Zg->i[cnt] = W->i[k];
      Zg->x[cnt] = W->x[k];
      cnt++;
    } 
    for(k = 0; k < qi; k++) Zg->p[k] = W->p[si+k] - W->p[si];
    Zg->p[qi] = cnt;

    // Meyer 1997: eqn. 20 [also f(theta) in Johnson & Thompson 1995 eqn 9c)]
    // Knight 2008: eqn 2.32-2.35 and 3.11-3.13 (p.34)
    //TODO for covariances see Johnson & Thompson 1995 eqn 11b
    //// B[,g] = Z_g %*% sln[si:ei, 1] %*% Ginv[g] 
    ////// See difference in order between Johnson & Thompson (e.g., eqn. 11b) and Meyer 1997 eqn 20
    ////// Actually calculates t(sln) %*% Z_g for computational ease (go by columns)
    //////// Also include multiplication with Ginv
    ////////// because matrix multiplication identity holds:
    ////////// (Ginv * (r_1 * c_1 + r_2 * c_2) == Ginv * r_1 * c_1 + Ginv * r_2 * c_2)

    ////// alternatively cs_gaxpy(Zg, sln_g, Btmp);  /* y = A*x+y */
    for(i = 0; i < qi; i++){
      for(k = Zg->p[i]; k < Zg->p[i+1]; k++){
        B->x[B->p[g] + Zg->i[k]] += Zg->x[k] * BLUXs->x[si + i] * Ginv[g]->x[0];  //TODO what to do if Ginv is a matrix?
      }
    }

    cs_spfree(Zg);
    si = ei + 1;
  }  // end for g

  //FIXME TODO Check what to do if more than 1 residual variance parameter
  for(k = 0; k < n; k++){
    B->x[B->p[nG] + k] += rory[k] * Rinv->x[0];  //TODO what to do if Rinv is a matrix?
  }
  tB = cs_transpose(B, 1);

  
  // Set up modified MME like the MMA of Meyer 1997 eqn. 23
  //// Substitute `B` instead of `y`

  // lambda: (co)variance ratios
  if(lambda == 1){
    BRHS = cs_multiply(tW, B);
    tBRinvB = cs_multiply(tB, B);

  //// NOT lambda: Straight (co)variances
  }else{
    BRHS = cs_multiply(tWKRinv, B);
    tBKRinv = cs_multiply(tB, KRinv);
    // next is actually `tBKRinvB`, want 1 name for this and when lambda=TRUE
    tBRinvB = cs_multiply(tBKRinv, B);
      cs_spfree(tBKRinv);
  }    // end if lambda=FALSE
  // double-transpose to get BRHS in correct order
  tBKRinv = cs_transpose(BRHS, 1);
    cs_spfree(BRHS);
  BRHS = cs_transpose(tBKRinv, 1);
    cs_spfree(tBKRinv);
  BRHSx = BRHS->x;



  // tBPB
  //// Johnson & Thompson 1995 eqns 8,9b,9c
  ////// (accomplishes same thing as Meyer 1997 eqns 22-23 for a Cholesky
  ////// factorization as Boldman & Van Vleck eqn 6 applied to AI calculation)
  // create temporary column of S as a 1 column matrix
  S = cs_spalloc(BRHS->m, BRHS->n, BRHS->m * BRHS->n, true, false);
    Sp = S->p;
    Sp[0] = 0;
  cnt = 0;
  for(k = 0; k < p; k++){
    g = 0;
    for(i = 0; i < BRHS->m; i++){
      // permute when make p_sln so skip `cs_ipvec` step below
      if(i == BRHS->i[BRHS->p[k] + g]){
        p_sln[Pinv[i]] = BRHSx[BRHS->p[k] + g];
        g++;
      } else{ 
        p_sln[Pinv[i]] = 0.0; 
      }
      Scol[i] = 0.0;
    }  // end for i
    
    cs_lsolve(Lc, p_sln);                         // x = L\x 
    cs_ltsolve(Lc, p_sln);		          // x = L'\x 
    cs_pvec(Pinv, p_sln, Scol, BRHS->m);          // b = P'*x 
  
    // put Scol into kth column of S
    for(i = 0; i < BRHS->m; i++){
      S->i[cnt] = i;
      S->x[cnt] = Scol[i];
      cnt++;
    }  // end for i
    Sp[k + 1] = cnt;

  }  // end for k in columns of S (ncol(S)=p)
  S->p = Sp;
  tS = cs_transpose(S, 1);


  // AI = -0.5 * (tBRinvB - tS %*% BRHS)
  //// cs_add() below = `alpha*A + beta*B` gives 1*tBRinvB + -1*tSBRHS
  tSBRHS = cs_multiply(tS, BRHS);
  AI = cs_add(tBRinvB, tSBRHS, 1.0, -1.0); 
  for(i = 0; i < AI->p[AI->n]; i++) AI->x[i] *= 0.5;
  if(lambda == 1) for(i = 0; i < AI->p[AI->n]; i++) AI->x[i] /= sigma2e;
  //////////////////////////////////////////////////////////////////////////////

  cs_spfree(tSBRHS);
  cs_spfree(tBRinvB);  
  cs_spfree(tS);
  cs_spfree(BRHS);
  cs_spfree(S);
  cs_spfree(tB);
  cs_spfree(B);
  cs_spfree(Rinv);

  delete [] p_sln;
  delete [] Scol;

 return (AI);
}



