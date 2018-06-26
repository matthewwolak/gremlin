#include "gremlincc.h"
/*******************************************************************/
/* 2 clock functions from SuiteSparse 5.1.0 by Tim Davis           */
static double tic (void) { return (clock () / (double) CLOCKS_PER_SEC) ; }
static double toc (double t) { double s = tic () ; return (CS_MAX (0, s-t)) ; }
/*******************************************************************/





extern "C"{  
void ugremlin(
	double *y,		// response
	int *ny,		// No. responses
	double *D,		// crossprod(y) <--> corner of M
	int *ndgeninv,		// Non-diagonal ginverse (geninv) matrices
//	int *nndgeninv,		// No. of non-diagonal geninv matrices
	int *dimZgs,		// dimensions of Zg matrices
	int *dimXZWG,		// dimensions of X, Z, W, & geninv matrices
//FIXME: if don't pass in X: rename and change indexing of this below
	int *nnzXWG,		// No. non-zeros in X, W, and geninvs
	int *iX,		// X
	int *pX,
	double *xX,
	int *iW,		// W
	int *pW,
	double *xW,
	int *igeninv,		// geninv (each geninv)
	int *pgeninv,
	double *xgeninv,
	double *rfxlL,		// Random effects contribution to log-Likelihood
	int *p,			// No. theta parameters
	int *nGR,		// No. G and R theta parameters
	int *dimGRs,		// dimensions of G and R matrices
	int *iGRs,		// GRs (G and R matrices of theta parameters)
	int *pGRs,		
	int *nnzGRs,		// No. non-zeroes in GRs	
	double *theta,		// thetas
	int *nnzBpinv,		// inverse Fixed effects prior matrix (no prior=0s)
	int *iBpinv,
	int *pBpinv,
	double *xBpinv,
	double *dLdtheta,	// gradient vector (1st deriv. of lL / thetas)
	double *lHvec,		// vector for lower triangle of information matrix (2nd deriv. of lL / thetas)
	double *sln,		// solution vector
        double *Cinv_ii,	// diagonal of Cinv: sampling variances of BLUXs
	double *r,		// residual vector
	double *itMat,		// parameter information at each iteration
	int *algit,		// algorithm for optimization at each iteration
	int *maxit,		// maximum No. iterations
	double *cctol,		// convergence criteria tolerances
	double *ezero,		// effective zero
	int *v,			// verbosity level of  output messages
	int *vit		// at what iterations to output messages
		
){

  cs	 *Bpinv, *W, *tW, *tWW, *tWWtmp, *cRinv, *tcRinv,
	 *Ctmp, *C, *tC, *Lcij, *Lc, *pCinv, *Cinv,
	 *RHS, *tCRHS, *CRHS, *RHSD, *M, *tmpBLUXs, *BLUXs, *R;
  css    *sLm;
  csn    *Lm;
  csi    Cn, *P;
  csi    *Pinv = new csi[dimXZWG[5]];

  int    nG = nGR[0], nR = nGR[1];

  cs*    *geninv = new cs*[nG];
  cs*    *G = new cs*[nG];
  cs*    *GcRinv = new cs*[nG];
  cs*    *GRinv = new cs*[nG];
  cs*    *invGRinv = new cs*[nG];
  cs* 	 *KGRinv = new cs*[nG];

  double t, T, took, dsLc, tyPy, logDetC, sigma2e, loglik, d, cc2, cc2d;

  int 	 g, i, k, rw, si, si2, ei,
	 itc = 0,
         dimM,
	 nminffx = ny[0] - dimXZWG[1],
	 nr = dimXZWG[3],
	 nminfrfx;

  int	 *cc = new int[5];

  char   tookUnits;





if(v[0] > 3) t = tic();


  d = 0.0;    // temporary for calculating change in theta (cc2) and EM residual

  //FIXME Do I need X? Delete if not
  // setup X: fixed-effects design matrix
//  X = cs_spalloc(dimXZWG[0], dimXZWG[1], nnzXWG[0], true, false);
//    for(k = 0; k < nnzXWG[0]; k++){
//      X->i[k] = iX[k];
//      X->x[k] = xX[k];
//    }
//    for(k = 0; k <= dimXZWG[1]; k++){
//      X->p[k] = pX[k];
//    }

  // Setup all ginverse (geninv) matrices
  //// (leave elements for diagonal/I matrices alone)
  si = 0; si2 = 0;
  for(g = 0; g < nG; g++){
    if(ndgeninv[g] == 1){
      dimM = dimXZWG[6+g];
      geninv[g] = cs_spalloc(dimM, dimM, nnzXWG[2+g], true, false);
      for(k = 0; k < nnzXWG[2+g]; k++){
        geninv[g]->i[k] = igeninv[si+k];
        geninv[g]->x[k] = xgeninv[si+k];
      }
      for(k = 0; k <= dimM; k++){
        geninv[g]->p[k] = pgeninv[si2+k];
      }
      si += nnzXWG[2+g]; si2 += dimM+1;
    }  // End if non-diagonal
  }  // End for g


  // Setup inverse matrix of Fixed effects prior (Bpinv)
  //// See Schaeffer 1991's summary of Henderson's development on this
  //// empty/zeroes if no prior specified
  Bpinv = cs_spalloc(dimXZWG[1], dimXZWG[1], nnzBpinv[0], true, false);
  for(k = 0; k < nnzBpinv[0]; k++){
    Bpinv->i[k] = iBpinv[k];
    Bpinv->x[k] = xBpinv[k];
  }
  for(k = 0; k <= dimXZWG[1]; k++){
    Bpinv->p[k] = pBpinv[k];
  }  
  //TODO: do I need to cs_sprealloc() to get rid of zeroes explicitly stored?


  // 1 Create mixed model array (M) and coefficient matrix of MME (C)
  //// quadratic at bottom-right (Meyer & Smith 1996, eqn 6)
  //// setup W=[X Z]
  W = cs_spalloc(dimXZWG[4], dimXZWG[5], nnzXWG[1], true, false);
    for(k = 0; k < nnzXWG[1]; k++){
      W->i[k] = iW[k];
      W->x[k] = xW[k];
    }
    for(k = 0; k <= dimXZWG[5]; k++){
      W->p[k] = pW[k];
    }
  tW = cs_transpose(W, true);
  tWW = cs_multiply(tW, W); //TODO add statement to include `Rinv`
    // Now take transpose of transpose to correctly order (don't ask why)
    tWWtmp = cs_transpose(tWW, true);
    cs_spfree(tWW);
    tWW = cs_transpose(tWWtmp, true);
    cs_spfree(tWWtmp);

  // setup 0 initialized RHS 1-column matrix
  //// ?Potentially: having explicit 0s in RHS will ensure M ordered with RHS last
  ////TODO check this last statement
  //// XXX KEEP explicit 0s for BLUXs cholesky/triangular solve to work
  RHS = cs_spalloc(dimXZWG[5], 1, dimXZWG[5], true, false);
  tmpBLUXs = cs_spalloc(dimXZWG[5], 1, dimXZWG[5], true, false);  // *sln matrix
  BLUXs = cs_spalloc(dimXZWG[5], 1, dimXZWG[5], true, false);    // sln matrix
  for(k = 0; k < dimXZWG[5]; k++){
    RHS->i[k] = k;
    tmpBLUXs->i[k] = k;
    BLUXs->i[k] = k;
    RHS->x[k] = 0.0;
    tmpBLUXs->x[k] = 0.0;
    BLUXs->x[k] = 0.0;
  }
  RHS->p[0] = 0; RHS->p[1] = dimXZWG[5];
  tmpBLUXs->p[0] = 0; tmpBLUXs->p[1] = dimXZWG[5];
  BLUXs->p[0] = 0; BLUXs->p[1] = dimXZWG[5];


  // setup matrix for residuals
  R = cs_spalloc(dimXZWG[5], 1, dimXZWG[5], true, false);
  for(k = 0; k < ny[0]; k++){
    R->i[k] = k;
  }
  R->p[0] = 0; R->p[1] = ny[0];





  // setup R matrix
  //// Assumes, just 1 R matrix 
  si = 0; for(g = 0; g < nG; g++) si += nnzGRs[g];
  dimM = dimGRs[nG];
  cRinv = cs_spalloc(dimM, dimM, nnzGRs[nG], true, false);
  for(k = 0; k < nnzGRs[nG]; k++){
    cRinv->i[k] = iGRs[si+k];
    cRinv->x[k] = 1.0 / sqrt(theta[si+k]);
  }
  for(k = 0; k <= dimM; k++){
    cRinv->p[k] = k*dimM;
  }
  //TODO Transform starting parameters to 'nu' scale
  //// cholesky of covariance matrices, then log of transformed diagonals
  //////TODO reference for above method (Meyer?)

/*TODO re-instate after make multivariate (also need an R matrix setup just like cRinv is above
  if(cR == NULL){
    error("R-structure %i starting value(s) is/are improperly specified: check that all eigenvalues (`eigen(R)$values`) > 0 and that the cholesky decomposition can be formed (`chol(R)`)\n", nR);
  }
*/
  tcRinv = cs_transpose(cRinv, true);






  // setup G matrices
  //// FIXME: with strategy of factoring out univariate sigma2E (see Meyer 1989)
  si = 0;
  for(g = 0; g < nG; g++){
    dimM = dimGRs[g];
    G[g] = cs_spalloc(dimM, dimM, nnzGRs[g], true, false);
    for (k = 0; k < nnzGRs[g]; k++){        
      G[g]->i[k] = iGRs[si+k];
      G[g]->x[k] = theta[si+k];
    }
    si += nnzGRs[g];
    for(k = 0; k<= dimM; k++){
      G[g]->p[k] = k*dimM;
    }
    GcRinv[g] = cs_multiply(tcRinv, G[g]);  // Next two lines from Meyer 1991, p.77
    GRinv[g] = cs_multiply(GcRinv[g], cRinv);
    invGRinv[g] = cs_inv(GRinv[g]);
  }





  //////////////////////////////////////////////////////////////////////////////
  // 1c Now make coefficient matrix of MME
  //// invert G_i and form Kronecker product with ginverse element (i.e., I or geninv)
  for(g = 0; g < nG; g++){
    if(ndgeninv[g] == 0){   // Diagonal matrix kronecker product: invGRinv %x% I
      KGRinv[g] = cs_kroneckerI(invGRinv[g], dimXZWG[6+g]);
    }else{  // generalized inverse kronecker product: invGRinv %x% geninv
      KGRinv[g] = cs_kroneckerA(invGRinv[g], geninv[g]);
    }
  }

  //// Construct C
  if(nG > 0){
    Ctmp = cs_omega(KGRinv, nG, Bpinv); 
    C = cs_add(tWW, Ctmp, 1.0, 1.0);
  }else{
    C = cs_add(tWW, Bpinv, 1.0, 1.0);
  }
  Cn = C->n;




if(v[0] > 3){
  took = toc(t); 
  Rprintf("%6.4f sec. (CPU clock): initial cpp setup (to get C)\n", took);
  t = tic();
}




  // Create Mixed Model Array (M) matrix
  cs_gaxpy(tW, y, RHS->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x`
  RHSD = cs_spalloc(dimXZWG[5]+1, 1, RHS->nzmax+1, true, false);
  si = 0;
  for(k = 0; k < RHS->nzmax; k++){
    RHSD->i[k] = RHS->i[k];
    RHSD->x[k] = RHS->x[k];
    si++;
  }
  RHSD->i[si] = dimXZWG[5];
  RHSD->x[si] = D[0];
  RHSD->p[0] = 0; RHSD->p[1] = si + 1;
 
  tCRHS = cs_cbind(C, RHS);
  CRHS = cs_transpose(tCRHS, true);
  cs_spfree(tCRHS);
  M = cs_cbind(CRHS, RHSD);
  cs_spfree(CRHS); cs_spfree(RHSD);

  //1d Find best order of MMA/C/pivots!
  //// Graser et al. 1987 (p1363) when singular C, |C| not invariant to order
  //// of C, therefore same order (of pivots) must be used in each loglik iteration
  //// Hadfield (2010, MCMCglmm paper App B): details on chol, ordering, and updating
  
  // Symbolic Cholesky factorization of M
  sLm = cs_schol(1, M);
  if(sLm == NULL){
      error("FAILED: symbolic Cholesky factorization of Mixed Model Array (`M`)\n");
  }
  if(sLm->pinv[Cn] != Cn){
    Rprintf("sLm->pinv[%i]=%i\n", Cn, sLm->pinv[Cn]);
    error("Ordering of M has not left RHS as last row/column\n");
  }
  // Allocate permutation matrix for C (all but last entry)
  for(k = 0; k < Cn; k++) Pinv[k] = sLm->pinv[k];
  P = cs_pinv(Pinv, Cn);


if(v[0] > 3){
  took = toc(t);
  Rprintf("  %6.4f sec. (CPU clock): initial cpp cs_schol(M)\n", took);
  t = tic();
}





  // Create variables needed for log-likelihood calculations
  // Used for log(|R|) and log(|G|) <-- Meyer 1989 (univar.) & 1991 (multivar.)
  int	 *rfxlvls = new int[nG];
    for(g = 0; g < nG; g++){
      rfxlvls[g] = dimZgs[2*g+1];
    }
    nminfrfx = nminffx - nr; 



if(v[0] > 3){
  took = toc(t);
  Rprintf("  %6.4f sec. (CPU clock): rest of initial cpp setup\n", took);
}

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //XXX									     XXX
  /*		REML ITERATIONS			REML ITERATIONS		      */
  //XXX									     XXX
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  for(i = 0; i < maxit[0]; i++){
    T = tic();
    if(v[0] > 0 && i%vit[0] == 0){
      Rprintf("  %i of max %i\n", i+1, maxit[0]); //TODO TIME of DAY?
    }

    if(i > 0){
      cs_spfree(tcRinv);
      for(g = 0; g < nG; g++){
        cs_spfree(GcRinv[g]); cs_spfree(GRinv[g]); cs_spfree(invGRinv[g]);
      }
      cs_spfree(C);
      cs_nfree(Lm);



      // setup R matrix
      //// Assumes, just 1 R matrix 
      si = 0; for(g = 0; g < nG; g++) si += nnzGRs[g];
      dimM = dimGRs[nG];
      for(k = 0; k < nnzGRs[nG]; k++){
        cRinv->x[k] = 1.0 / sqrt(theta[si+k]);
      }

    //TODO Transform starting parameters to 'nu' scale (don't need for EM algorithm)
    //// cholesky of covariance matrices, then log of transformed diagonals
    //////TODO reference for above method (Meyer?)

/*TODO re-instate after make multivariate
  if(cR == NULL){
    error("R-structure %i starting value(s) is/are improperly specified: check that all eigenvalues (`eigen(R)$values`) > 0 and that the cholesky decomposition can be formed (`chol(R)`)\n", nR);
  }
*/
      tcRinv = cs_transpose(cRinv, true);





      // setup G matrices
      //// FIXME: with strategy of factoring out univariate sigma2E (see Meyer 1989)
      si = 0;
      for(g = 0; g < nG; g++){
        dimM = dimGRs[g];
        for (k = 0; k < nnzGRs[g]; k++){        
          G[g]->x[k] = theta[si+k];
        }
        si += nnzGRs[g];
        GcRinv[g] = cs_multiply(tcRinv, G[g]);  // Next two lines from Meyer 1991, p.77
        GRinv[g] = cs_multiply(GcRinv[g], cRinv);
        invGRinv[g] = cs_inv(GRinv[g]);
      }





      //////////////////////////////////////////////////////////////////////////
      // 1c Now make coefficient matrix of MME
      //// invert G_i and form Kronecker product with ginverse element (i.e., I or geninv)
      for(g = 0; g < nG; g++){
        if(ndgeninv[g] == 0){   // Diagonal matrix kronecker product: invGRinv %x% I
          cs_kroneckerIupdate(invGRinv[g], dimXZWG[6+g], KGRinv[g]);
        }else{  // generalized inverse kronecker product: invGRinv %x% geninv
          cs_kroneckerAupdate(invGRinv[g], geninv[g], KGRinv[g]);
        }
      }

      //// Construct C
      if(nG > 0){
        cs_omegaupdate(KGRinv, nG, Bpinv, Ctmp); 
        C = cs_add(tWW, Ctmp, 1.0, 1.0);
      }else{
        C = cs_add(tWW, Bpinv, 1.0, 1.0);
      }
      // order rows within columns (necessary for updating M)
      tC = cs_transpose(C, true);
      cs_spfree(C);
      C = cs_transpose(tC, true);
      cs_spfree(tC);

      // update M
      // non-zero pattern hasn't changed - nor has last row & column
      for(k = 0; k < Cn; k++){
        for(g = C->p[k], rw = M->p[k]; g < C->p[k+1]; g++, rw++){
          M->x[rw] = C->x[g];
        }
      }


    }  // End if i>0



if(v[0] > 3){
  took = toc(t);
  Rprintf("  %6.4f sec. (CPU clock): cpp REML i=%i setup\n", took, i);
  t = tic();
}




    Lm = cs_chol(M, sLm);  //TODO update if i>0? (cs_updown)
//    Lm = cs_chol(M, sLm);  //TODO update if i>0? (cs_updown)
    if(Lm == NULL){
      error("Mixed Model Array singular: possibly caused by a bad combination of G and R (co)variance parameters\n");
    }


if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec. (CPU clock): cpp REML i=%i cs_chol(M)\n", took, i);
  t = tic();
}



    // Obtain Cholesky factor of C (Lc) from Lm (first Cn rows & columns of Lm)
    //// Assume, same non-zero pattern each iteration (i.e., order or rows same)
    if(i == 0){
      // form a triplet matrix
      //// allocate max. entries according to non-zeroes in m-1 columns of M 
      //// does include the non-zeroes from last row of M
      Lcij = cs_spalloc(Cn, Cn, Lm->L->p[Cn], true, true);
      k = 0; rw = 0;
      for(g = 0; g < Lm->L->p[Cn]; g++){
        if(g >= Lm->L->p[k+1]) k++;
        if(Lm->L->i[g] != Cn){
          Lcij->i[rw] = Lm->L->i[g];
          Lcij->p[rw] = k;
          Lcij->x[rw] = Lm->L->x[g];
          rw++;
        }
      Lcij->nz = rw;
      }

    } else{
      rw = 0;
      for(g = 0; g < Lm->L->p[Cn]; g++){
        if(Lm->L->i[g] != Cn){
          Lcij->x[rw] = Lm->L->x[g];
          rw++;
        }
      }
    }    
    // Convert triplet matrix into sparse compressed column matrix
    Lc = cs_compress(Lcij);


if(v[0] > 3){
  took = toc(t);
  Rprintf("\t    %6.4f sec. (CPU clock): cpp REML i=%i chol(C) from chol(M)\n", took, i);
  t = tic();
}







  



    /*		XXX		XXX		XXX		XXX	*/
    pCinv = cs_chol2inv(Lc);
    if(pCinv == NULL){
      error("Error inverting C\n");
    }

    Cinv = cs_permute(pCinv, P, Pinv, 1);  //<-- permute values 0/1=FALSE/TRUE
    // Do t(t(Cinv)) to order correctly
    cs_spfree(pCinv);   
    pCinv = cs_transpose(Cinv, true);
    cs_spfree(Cinv);
    Cinv = cs_transpose(pCinv, true);
    cs_spfree(pCinv);


if(v[0] > 3){
  took = toc(t);
  Rprintf("\t    %6.4f sec. (CPU clock): cpp REML i=%i chol2inv(Lc) to Cinv\n", took, i);
  t = tic();
}

    ////////////////////////////////////////////////////////////////////////////
    // 5 record log-like, check convergence, & determine next varcomps to evaluate  
    //// 5a determine log(|C|) and y'Py
    ////// Meyer & Smith 1996, eqns 12-14 (and 9)
    //// Also see Meyer & Kirkpatrick 2005 GSE. eqn. 18: if cholesky of MMA = LL'
    dsLc = 0.0; tyPy = 0.0; loglik = 0.0;  // set diagonal sum of Lc, tyPy, and loglik back to 0
    // Meyer & Smith 1996, eqn. 13
    for(k = 0; k < Cn; k++){
      rw = Lm->L->p[k];
      dsLc += log(Lm->L->x[rw]);
    }
    logDetC = 2.0 * dsLc;
    // Meyer & Smith 1996, eqn. 14
    tyPy += Lm->L->x[Lm->L->p[Cn]];
      tyPy *= tyPy;

    sigma2e = tyPy / nminffx;             // residual variance
    loglik += nminffx + logDetC;         // nminffx is tyPy/sigma2e, simplified b/c sigma2e=tyPy/nminffx 
    
    // log(|R|) contribution TODO Assumes X of full rank
    //// Alternatively: if Rinv NOT factored out of MMA `loglik += ny*log(start$R)`
    loglik += nminfrfx * log(sigma2e);
    // log(|G|) contribution
    //// FIXME only works for independent random effects
    //// See also alternative (using starting value for residual when factored Rinv)
    if(nG > 0){
      for(g = 0; g < nG; g++){
        //FIXME below assumes only one entry in GRinv[g]
        loglik += rfxlvls[g] * log(GRinv[g]->x[0] * sigma2e);
      }
    }
    loglik += rfxlL[0];
    loglik *= -0.5;



if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec. (CPU clock): cpp REML i=%i log-likelihood calc.\n", took, i);
  t = tic();
}

    //XXX	**** END LOG-LIKELIHOOD CALCULATION **** 		XXX
    //XXX	**** END LOG-LIKELIHOOD CALCULATION **** 		XXX
    /**************************************************************************/




    ////////////////////////////////////////////////////////////////////////////
    // solve MME for BLUEs/BLUPs
////TODO Might not need to solve for BLUEs (see Knight 2008 for just getting BLUPs eqn 2.13-15
    //// see Mrode 2005 chapter
    //XXX tmpBLUXs must have EXPLICIT 0s, else triangular solve = incorrect answer
    cs_ipvec (Pinv, RHS->x, tmpBLUXs->x, Cn) ;   /* x = P*b */
    cs_lsolve (Lc, tmpBLUXs->x) ;                /* x = L\x */
    cs_ltsolve (Lc, tmpBLUXs->x) ;               /* x = L'\x */
    cs_pvec (Pinv, tmpBLUXs->x, BLUXs->x, Cn) ;  /* b = P'*x */

if(v[0] > 3){
  took = toc(t);
  Rprintf("\t    %6.4f sec. (CPU clock): cpp REML i=%i sln forward/back solve with chol(C)\n", took, i);
  t = tic();
}


    // calculate residuals as r = y - W %*% sln
    // put y in R->x to be over-written
    for(k = 0; k < ny[0]; k++){
      R->x[k] = -1.0 * y[k];
    }
    //// below rearranges to -r = W %*% sln - y
    cs_gaxpy(W, BLUXs->x, R->x);   //  y = A*x+y
    for(k = 0; k < ny[0]; k++){
      R->x[k] *= -1.0;
    }



if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec. (CPU clock): cpp REML i=%i sln/r calc.\n", took, i);
  t = tic();
}



    for(k = 0; k < p[0]; k++){
      itMat[itc] += theta[k];
      itc++;
    }
    itMat[itc] += sigma2e; itc++;
    itMat[itc] += tyPy;    itc++;
    itMat[itc] += logDetC; itc++;
    itMat[itc] += loglik;  itc++;





if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec. (CPU clock): cpp REML i=%i itMat recording took\n", took, i);
  t = tic();
}


    ////////////////////////////////////////////////////////////////////////////
    // 5c check convergence criteria
    //// Knight ('08 ch. 6): Searle et al.'92 & Longford '93 discuss diff. crit. types
    //// Appendix 2 of WOMBAT help manual for 4 criteria specified
    for(k = 0; k <= 4; k++) cc[k] = 0;     // Keep track of sum in last position
    d = 0.0; cc2 = 0.0; cc2d = 0.0;
    if(i > 0){
      // Change in log-likelihood
      cc[0] += (itMat[itc-6-p[0]] - loglik) < cctol[0];
        cc[4] += cc[0];
      // Change in theta: wombat eqn. A.1 also Knight 2008 eqn. 6.1
      for(k = 0; k < p[0]; k++){
        d = theta[k] - itMat[itc-(9+2*p[0])+k];
        cc2 += d*d;
        cc2d += theta[k] * theta[k];
      }
      cc[1] += sqrt(cc2 / cc2d) < cctol[1];
        cc[4] += cc[1];
      // 3 & 4 are only for optimzation algorithms which produce derivatives
      if(algit[i] > 0){    // 0=EM, 1=AI
        // Norm of gradient vector: wombt eqn. A.2
          //TODO Does this go here or maybe after AI
          //// Does this step happen for last AI matrix (i-1) or current (i)?
          //  cc[4] += cc[2];
        // Newton decrement: wombat eqn A.3 and Boyd & Vandenberghe 2004
          //TODO
          //  cc[4] += cc[3];
      }
    }  // end convergence criteria if i>0




if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec. (CPU clock): cpp REML i=%i convergence crit. calc.\n", took, i);
  t = tic();
}



    ////////////////////////////////////////////////////////////////////////////
    // 5d Determine next (co)variance parameters to evaluate: REML NOT CONVERGED
    // Meyer and Smith 1996 for algorithm using derivatives of loglik
    //// eqn 15-18 (+ eqn 33-42ish) for derivatives of tyPy and logDetC
    //// Smith 1995 for very technical details
//XXX trace of a matrix corresponds to the derivative of the determinant
    if(cc[4] < 2){
      // Expectation Maximization
      // EM refs: Hofer 1998 eqn 10-12
      // XXX note Hofer eqn 12 missing sigma2e in last term of non-residual formula
      ////XXX see instead Mrode 2005 (p. 241-245)
      if(algit[i] == 0){
        if(v[0] > 1 && i%vit[0] == 0) Rprintf("\t EM to find next theta\n");
        if(!cs_em(BLUXs, theta, nG, rfxlvls, dimXZWG[1], ndgeninv, geninv, Cinv)){
          error("Unusccessful EM algorithm in iteration %i\n");
        }
        // Calculate EM for residual:
        //// crossprod(y, r) / nminffx
        d = 0.0;
        for(k = 0; k < ny[0]; k++) d += y[k] * R->x[k]; // assumes R has explicit 0s
        theta[nG] = d / nminffx;
      }  // end EM

      // Average Information
      // Appendix 5 WOMBAT manual:how to modify AI matrix ensure improvements of logLik
      if(algit[i] == 1){
        if(v[0] > 1 && i%vit[0] == 0) Rprintf("\t AI to find next theta\n");
        //TODO do I need to check convergence criteria here (i.e., cc[3:4])
      }  // end AI


    }  // end if REML did not converge







    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    took = toc(T);                 // Capture cpu clock time for i REML iteration
    // V=1 LEVEL of OUTPUT
    if(v[0] > 0 && i%vit[0] == 0){ 
      Rprintf("\t\tlL:%6.6f", loglik);
      Rprintf("\ttook %6.3f sec. (CPU clock)\n", took); //TODO format units if >60 (and do for all Rprintf(took))
      // To format units see const *char in: http://www.cplusplus.com/reference/cmath/round/

      // V=2 LEVEL of OUTPUT
      if(v[0] > 1){
        Rprintf("\n");
        // output names of variables
        for(g = 0; g < p[0]; g++){
          Rprintf("\tV%i", g+1);
        }
        Rprintf("\tsigma2e\ttyPy\tlogDetC\n");
        // output underlines"
        for(g = 0; g < p[0]; g++){
          if(g < 10) Rprintf("\t--"); else Rprintf("\t---");
        }
        Rprintf("\t------\t----\t--------\n");
        // output variable values
        for(g = 0; g < p[0]; g++){
          Rprintf("\t%6.4f", itMat[itc-(4+p[0])+g]);
        }
        Rprintf("\t%6.4f\t%6.4f\t%6.4f\n", sigma2e, tyPy, logDetC);

        // output convergence criteria
        Rprintf("\tConvergence crit:");
          for(k = 0; k < 4; k++) Rprintf("%4i", cc[k]);
          Rprintf("\n");
        
        // V=3 LEVEL of OUTPUT
        if(v[0] > 2){
          Rprintf("\n\n"); //TODO
        }  // end if v>2
      }  // end if v>1
    }  // end if v>0


    itMat[itc] += round(took*10) / 10;              // gives 1 decimal place
    itc++;
 



   
    // Determine if model has converged
    //FIXME: change number of parameters that must be true as add criteria
    if(cc[4] > 1){ // FIXME: change number of parameters that must be true
      Rprintf("\n\nREML converged\n\n");      
      break;
    }


  }  // end i for loop 

  //////////////////////////////////////////////////////////////////////////////
  //XXX									     XXX
  /*		END reml iterations		END reml iterations	      */
  //XXX									     XXX
  //////////////////////////////////////////////////////////////////////////////



if(v[0] > 3) t = tic();

  // Pass information to R
  //// Solution vector
  for(k = 0; k < dimXZWG[5]; k++) sln[k] += BLUXs->x[k];
  //// diagonals of Cinv (sln sampling variances) 
  for(k = 0; k < Cn; k++){
    for(j = Cinv->p[k]; j < Cinv->p[k+1]; j++){
      if(Cinv->i[j] == k){
        Cinv_ii[k] += Cinv->x[j];
        break; 
      }  // end if
    }  // end for j
  }  // end for k
  //// Residual vector
  for(k = 0; k < ny[0]; k++) r[k] += R->x[k];
  //TODO
  //// return AI and dLdtheta

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //  cs_spfree(X);
  cs_spfree(Bpinv);
  cs_spfree(W); cs_spfree(tW); cs_spfree(tWW);
  cs_spfree(cRinv); cs_spfree(tcRinv);
  cs_spfree(Ctmp); cs_spfree(C); cs_spfree(Cinv);
  cs_spfree(Lcij); cs_spfree(Lc);
  cs_spfree(RHS); cs_spfree(tmpBLUXs); cs_spfree(BLUXs); cs_spfree(R);

  cs_sfree(sLm);

  cs_nfree(Lm);

//
  for(g = 0; g < nG; g++){
    if(ndgeninv[g] == 1) cs_spfree(geninv[g]);
    cs_spfree(G[g]); cs_spfree(GcRinv[g]); cs_spfree(GRinv[g]);
    cs_spfree(invGRinv[g]); cs_spfree(KGRinv[g]);
  } 
//
  delete [] geninv;
  delete [] G; delete [] GcRinv; delete [] GRinv;
  delete [] invGRinv; delete [] KGRinv;
//
  delete [] Pinv; delete [] P;
  delete [] rfxlvls;
  delete [] cc;


if(v[0] > 3){
  took = toc(t);
  Rprintf("%6.4f sec. (CPU clock): cpp post-REML freeing-up\n", took);
}

}
}

