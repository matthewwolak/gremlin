#include "gremlin.h"
#include "gremlincc.h"

/* returns log-likelihood if successful, 0 if not
   Replaces all objects that are updated based on changed `nu` values */
csn *reml(csi n, csi *dimZWG, csi nG, csi p, double *y,
	cs *Bpinv, cs *W, cs *tW, csi *rfxlvls, double rfxlL,
	cs *R, cs *Rinv, cs **G, cs **Ginv, csi *ndgeninv, cs **geninv,
	cs *KRinv, cs **KGinv, cs *tWKRinv, cs *tWKRinvW, cs *Ctmp,
	cs *RHS, cs *tmpBLUXs, cs *BLUXs, double *res,
	css *sLc,
	double *tyPy, double *logDetC, double *sigma2e,
	double tyRinvy, // lambda=TRUE same every iteration else 0.0 when FALSE
        int nminffx,  // lambda=TRUE else 0
	double *loglik,
	csi i, csi v, csi vitout, csi lmbda 
){

  int     g, k, rw, nminfrfx;
  double  t[2], took, dsLc;  
  cs	  *tyRinv, *C;
  csn     *Lc;
  
  if(!CS_CSC(W)) return (0);

  if(v > 3) simple_tic(t);
  *loglik = 0.0;
  
  nminfrfx = nminffx - dimZWG[1]; 
  

  //////////////////////////////////////////////////////////////////////////
  // 1 Setup to create coeficient matrix of MME (C)
  if(lmbda == 0){
    // setup tyRinv [t(y) %*% KRinv] to receive output from cs_gaxpy
    tyRinv = cs_spalloc(1, dimZWG[2], dimZWG[2], true, false);
      for(k = 0; k < dimZWG[2]; k++){
        tyRinv->i[k] = 0;
        tyRinv->p[k] = k;
        tyRinv->x[k] = 0.0;
      }
      tyRinv->p[dimZWG[2]] = dimZWG[2];
    // ASSUME cs_gaxpy does KRinv %*% y, which gives correct values for tyRinv@x
    cs_gaxpy(KRinv, y, tyRinv->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x` 
    for(k = 0; k < n; k++){
      tyRinvy += tyRinv->x[k] * y[k];  
    }
    
    // Next creates RHS
    //// Meyer '97 eqn 11
    //// (Note different order of RHS from Meyer '89 eqn 6; Meyer '91 eqn 4)
    // Reset RHS to 0.0 so cs_gaxpy works right
    if(i > 0) for(k = 0; k < dimZWG[3]; k++) RHS->x[k] = 0.0;
    cs_gaxpy(tWKRinv, y, RHS->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x` 
  }  // end if lambda=FALSE 
  

   
  // 1c Now make coefficient matrix of MME
  //// form Kronecker products for each G and ginverse element (i.e., I or geninv)
  ////// Update; originally made first iteration (i=0)
  if(i > 0){
    for(g = 0; g < nG; g++){
      if(ndgeninv[g] == 0){   // Diagonal matrix kronecker product: Ginv %x% I
        cs_kroneckerIupdate(Ginv[g], dimZWG[4+g], KGinv[g]);
      } else{  // generalized inverse kronecker product: Ginv %x% geninv
          cs_kroneckerAupdate(Ginv[g], geninv[g], KGinv[g]);
        }
    }
  }  // end if i>0
  
  //// Construct C
  if(nG > 0){
    cs_omegaupdate(KGinv, nG, Bpinv, Ctmp); 
    C = cs_add(tWKRinvW, Ctmp, 1.0, 1.0);
  } else{
      C = cs_add(tWKRinvW, Bpinv, 1.0, 1.0);
    }

if((i == 0) && (v > 3)){
  took = simple_toc(t); 
  Rprintf("  %6.4f sec.: initial cpp setup to get C\n", took);
  simple_tic(t);
}



  
  // Update Cholesky factorization of C
  Lc = cs_chol(C, sLc);  //TODO update when i>0? (cs_updown), but then initially make Lc when i=0
  if(Lc == NULL){
    Rprintf("\nreml: Mixed Model Coefficient matrix (C) singular, possibly caused by a bad combination of G and R (co)variance parameters");

    // Cholesky factorization did not work
    //// set log-likelihood to 0.0 to flag the error in gremlin.cc then exit
      *loglik = 0.0;
      cs_spfree(C);
      if(lmbda == 0) cs_spfree(tyRinv);

    return (Lc);
  }

if(v > 3){
  took = simple_toc(t);
  Rprintf("\n    %6.4f sec.: cpp REML i=%i cs_chol(C)", took, i);
  simple_tic(t);
}






  ////////////////////////////////////////////////////////////////////////////
  // solve MME for BLUEs/BLUPs
  //// Do now, because need solutions as part of log-likelihood calc. (for tyPy)
  //XXX tmpBLUXs must have EXPLICIT 0s, else triangular solve = incorrect answer
  cs_ipvec(sLc->pinv, RHS->x, tmpBLUXs->x, C->n) ;   /* x = P*b */
  cs_lsolve(Lc->L, tmpBLUXs->x) ;                    /* x = L\x */
  cs_ltsolve(Lc->L, tmpBLUXs->x) ;                   /* x = L'\x */
  cs_pvec(sLc->pinv, tmpBLUXs->x, BLUXs->x, C->n) ;  /* b = P'*x */

if(v > 3){
  took = simple_toc(t); 
  Rprintf("\n\t    %6.4f sec.: cpp REML i=%i sln forward/back solve with chol(C)", took, i);
  simple_tic(t);
}


  // calculate residuals as r = y - W %*% sln
  // put y in r to be over-written
  for(k = 0; k < n; k++){
    res[k] = -1.0 * y[k];
  }
  //// below rearranges to -r = W %*% sln - y
  cs_gaxpy(W, BLUXs->x, res);   //  y = A*x+y
  for(k = 0; k < n; k++){
    res[k] *= -1.0;
  }


if(v > 3){
  took = simple_toc(t);
  Rprintf("\n    %6.4f sec.: cpp REML i=%i sln/r calc.", took, i);
  simple_tic(t);
}



  ////////////////////////////////////////////////////////////////////////////
  // 5 record log-likelihood, check convergence, & determine next varcomps  
  //// 5a determine log(|C|) and y'Py
  // Boldman and Van Vleck 1991 eqn. 6 (tyPy) and 7 (log|C|)    
  // Meyer 1997 eqn. 13 (tyPy)

  // set diagonal sum of Lc (`L`), tyPy, and loglik back to 0/starting values
  *tyPy = tyRinvy;
  dsLc = 0.0; 
  // 'by hand' calculate `tyPy = tyRinvy -  crossprod(BLUXs, RHS)`
  for(k = 0; k < dimZWG[3]; k++){
    *tyPy -= BLUXs->x[k] * RHS->x[k];
  }
    
  // calculate logDetC from diagonals of Lc (`L`)
  for(k = 0; k < C->n; k++){
    rw = Lc->L->p[k];
    dsLc += log(Lc->L->x[rw]);
  }
  *logDetC = 2.0 * dsLc;

  // Factored out residual variance (only for lambda scale)
  if(lmbda == 1) *sigma2e = *tyPy / nminffx;

  // V=2 LEVEL of OUTPUT
  if(v > 2 && vitout == 0){
    Rprintf("\n\tsigma2e\ttyPy\tlogDetC\n");
    Rprintf("\t-------\t----\t--------\n");
    Rprintf("\t%6.4f\t%6.4f\t%6.4f\n", *sigma2e, *tyPy, *logDetC);
  }  // end v > 2



  // Construct the log-likelihood (Meyer 1997 eqn. 8)
  //// (firt put together as `-2 log-likelihood`)
  // `log(|R|)`
  if(lmbda == 1){
    *loglik += nminfrfx * log(*sigma2e);
  } else{
      *loglik += n * log(R->x[0]);  //FIXME won't work when R (co)variances
    }
    if(v > 3) Rprintf("\n\t log|R|=%6.4f", *loglik);

  // `log(|G|)`
  //// Meyer 1997 eqn. 9 for lambda=FALSE equation
  //// NOTE: setting `sigma2e=1.0` above (when lambda=FALSE) means 1 line below
  //// FIXME only works for independent random effects
  if(nG > 0){
    for(g = 0; g < nG; g++){
      //FIXME below assumes only one entry in a G[g]
      *loglik += rfxlvls[g] * log(G[g]->x[0] * *sigma2e);
    }
  }
    if(v > 3) Rprintf("\n\t log|G| added=%6.4f", *loglik);
  *loglik += rfxlL;    // rfxIncContrib2loglik
    if(v > 3) Rprintf("\n\t\t log|G| rfxlL(%6.4f) added=%6.4f", rfxlL, *loglik);

  // log(|C|) + tyPy
  //// if lambda=TRUE then nminffx=tyPy/sigma2e simplified below
  ////// (because sigma2e=tyPy/nminffx)
  if(lmbda == 1){
    *loglik += *logDetC + nminffx;
  } else{
      *loglik += *logDetC + *tyPy;
    }
    if(v > 3) Rprintf("\n\t log|C|+tyPy added=%6.4f", *loglik);
    
  // Multiply by -0.5 to calculate `loglik` from `-2loglik`
  *loglik *= -0.5;
    if(v > 3) Rprintf("\n\t multiplied by -0.5=%6.4f", *loglik);


if(v > 3){
  took = simple_toc(t);
  Rprintf("\n    %6.4f sec.: cpp REML i=%i log-likelihood calc.", took, i);
  simple_tic(t);
}

  /////////////////////////////////////
  cs_spfree(C);
  if(lmbda == 0) cs_spfree(tyRinv);

 return (Lc);
}

