#include "gremlincc.h"
/*******************************************************************/
/* 2 clock functions from SuiteSparse 5.1.0 by Tim Davis           */
static double tic (void) { return ((clock_t) clock () / (double) CLOCKS_PER_SEC) ; }
static double toc (double t) { double s = tic () ; return (CS_MAX (0, s-t)) ; }
/*******************************************************************/





extern "C"{  
void ugremlin(
/*1 */	double *y,		// response
/*2 */	int *ny,		// No. responses
/*3 */	int *nminffx,		// No. observations - No. Fxd Fx
/*4 */	int *ndgeninv,		// Non-diagonal ginverse (geninv) matrices
/*5 */	int *dimZgs,		// dimensions of Zg matrices
/*6 */	int *dimZWG,		// dimensions of Z, W, & geninv matrices
/*7 */	int *nnzWG,		// No. non-zeros in W and geninvs
/*8 */	int *iW,		// W
/*9 */	int *pW,
/*10*/	double *xW,
/*11*/	int *igeninv,		// geninv (each geninv)
/*12*/	int *pgeninv,
/*13*/	double *xgeninv,
/*14*/	double *rfxlL,		// Random effects contribution to log-Likelihood
/*15*/	int *lambda,		// TRUE/FALSE lambda (variance ratios?)
/*16*/	int *p,			// No. nu parameters
/*17*/	int *nGR,		// No. G and R nu parameters
/*18*/	int *dimGRs,		// dimensions of G and R matrices
/*19*/	int *iGRs,		// GRs (G and R matrices of nu parameters)
/*20*/	int *pGRs,		
/*21 */	int *nnzGRs,		// No. non-zeroes in GRs	
/*22*/	double *nu,		// nus
/*23*/	int *nnzBpinv,		// inverse Fixed effects prior matrix (no prior=0s)
/*24*/	int *iBpinv,
/*25*/	int *pBpinv,
/*26*/	double *xBpinv,
/*27*/	double *dLdnu,		// gradient vector (1st deriv. of lL / nus)
/*28*/	double *AIvec,		// column-wise vector of average information matrix (2nd deriv. of lL / nu)
/*29*/	double *sln,		// solution vector
/*30*/  double *Cinv_ii,	// diagonal of Cinv: sampling variances of BLUXs
/*31*/	double *res,		// residual vector
/*32*/	double *itMat,		// parameter information at each iteration
/*33*/	int *algit,		// algorithm for optimization at each iteration
/*34*/	int *maxit,		// maximum No. iterations
/*35*/	double *cctol,		// convergence criteria tolerances
/*36*/	double *ezero,		// effective zero
/*37*/	int *v,			// verbosity level of  output messages
/*38*/	int *vit,		// at what iterations to output messages
/*39*/  int *sLcPinv		// empty Cholesky permutation vector
		
){

  cs	 *Bpinv, *W, *tW,
	 *R, *Rinv, *KRinv, *tyRinv, *tWKRinv, *tWKRinvW, *ttWKRinvW,
	 *Ctmp, *C,
	 *RHS, *tmpBLUXs, *BLUXs,
	 *AI, *AIinv;

  css    *sLc, *sLai;
  csn    *Lc, *Lai;

  int    nG = nGR[0];

  cs*    *geninv = new cs*[nG];
  cs*    *G = new cs*[nG];
  cs*    *Ginv = new cs*[nG];
  cs* 	 *KGinv = new cs*[nG];

  double t, T, took, dsLc, tyRinvy, tyPy, logDetC, sigma2e, loglik,
         d, cc2, cc2d;

  int 	 g, i, k, rw, si, si2, vitout,
	 itc = 0,
         aiformed = 0,
         dimM,      // GENERIC dimension of a matrix variable to be REUSED
	 nr = dimZWG[1],
	 nffx,
	 nminfrfx;

  int	 *rfxlvls = new int[nG];

  int	 *cc = new int[5];




if(v[0] > 3) t = tic();


  d = 0.0;    // temporary for calculating change in nu (cc2) and EM residual

  // Setup all ginverse (geninv) matrices
  //// (leave elements for diagonal/I matrices alone)
  si = 0; si2 = 0;
  for(g = 0; g < nG; g++){
    if(ndgeninv[g] == 1){
      dimM = dimZWG[4+g];	// dimensions for generalized inverses
      geninv[g] = cs_spalloc(dimM, dimM, nnzWG[1+g], true, false);
      for(k = 0; k < nnzWG[1+g]; k++){
        geninv[g]->i[k] = igeninv[si+k];
        geninv[g]->x[k] = xgeninv[si+k];
      }
      for(k = 0; k <= dimM; k++){
        geninv[g]->p[k] = pgeninv[si2+k];
      }
      si += nnzWG[1+g]; si2 += dimM+1;
    }  // End if non-diagonal
  }  // End for g


  nffx = ny[0] - nminffx[0];

  // Setup inverse matrix of Fixed effects prior (Bpinv)
  //// See Schaeffer 1991's summary of Henderson's development on this
  //// empty/zeroes if no prior specified
  Bpinv = cs_spalloc(nffx, nffx, nnzBpinv[0], true, false);
  for(k = 0; k < nnzBpinv[0]; k++){
    Bpinv->i[k] = iBpinv[k];
    Bpinv->x[k] = xBpinv[k];
  }
  for(k = 0; k <= nffx; k++){
    Bpinv->p[k] = pBpinv[k];
  }  
  //TODO: do I need to cs_sprealloc() to get rid of zeroes explicitly stored?



  // setup R matrix
  //// Assumes, just 1 R matrix 
  si = 0; for(g = 0; g < nG; g++) si += nnzGRs[g];
  dimM = dimGRs[nG];
  R = cs_spalloc(dimM, dimM, nnzGRs[nG], true, false);
  for(k = 0; k < nnzGRs[nG]; k++){
    R->i[k] = iGRs[si+k];
    R->x[k] = nu[si+k];
  }
  for(k = 0; k <= dimM; k++){
    R->p[k] = k*dimM;
  }

  Rinv = cs_inv(R);
/* 
could do a check to make sure R was inverted correctly:
  if(Rinv == NULL){
    error("R-structure %i starting value(s) is/are improperly specified: check that all eigenvalues (`eigen(R)$values`) > 0 and that the cholesky decomposition can be formed (`chol(R)`)\n", nR);
  }
*/





  // setup G matrices
  si = 0;
  for(g = 0; g < nG; g++){
    dimM = dimGRs[g];
    G[g] = cs_spalloc(dimM, dimM, nnzGRs[g], true, false);
    for (k = 0; k < nnzGRs[g]; k++){        
      G[g]->i[k] = iGRs[si+k];
      G[g]->x[k] = nu[si+k];
    }
    si += nnzGRs[g];
    for(k = 0; k<= dimM; k++){
      G[g]->p[k] = k*dimM;
    }
    Ginv[g] = cs_inv(G[g]);
  }






  //////////////////////////////////////////////////////////////////////////////
  // 1 Setup to create coefficient matrix of MME (C)
  // initialize W=[X Z]
  W = cs_spalloc(dimZWG[2], dimZWG[3], nnzWG[0], true, false);
    for(k = 0; k < nnzWG[0]; k++){
      W->i[k] = iW[k];
      W->x[k] = xW[k];
    }
    for(k = 0; k <= dimZWG[3]; k++){
      W->p[k] = pW[k];
    }
  tW = cs_transpose(W, true);





  // setup 0 initialized RHS 1-column matrix
  //// keep explicit 0s for cs_gaxpy to work in forming RHS
  RHS = cs_spalloc(dimZWG[3], 1, dimZWG[3], true, false);
    // since will be iterating over same dimensions, set up Solution matrices
    //// XXX KEEP explicit 0s for BLUXs cholesky/triangular solve to work
    tmpBLUXs = cs_spalloc(dimZWG[3], 1, dimZWG[3], true, false);  // *sln matrix
    BLUXs = cs_spalloc(dimZWG[3], 1, dimZWG[3], true, false);  // sln matrix
  for(k = 0; k < dimZWG[3]; k++){
    RHS->i[k] = k;
    tmpBLUXs->i[k] = k;
    BLUXs->i[k] = k;
    RHS->x[k] = 0.0;
    tmpBLUXs->x[k] = 0.0;
    BLUXs->x[k] = 0.0;
  }
  RHS->p[0] = 0; RHS->p[1] = dimZWG[3];
  tmpBLUXs->p[0] = 0; tmpBLUXs->p[1] = dimZWG[3];
  BLUXs->p[0] = 0; BLUXs->p[1] = dimZWG[3];



  // lambda: (co)variance ratios
  if(lambda[0] == 1){
    tWKRinvW = cs_multiply(tW, W);
    // Fill in RHS
    cs_gaxpy(tW, y, RHS->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x` 
    // form tyy (crossprod(y)), but call it tyRinvy for conciseness
    tyRinvy = 0.0;
    for(k = 0; k < ny[0]; k++){
      tyRinvy += y[k] * y[k];
    }

  //// NOT lambda: Straight (co)variances
  }else{
    // Rinv Kronecker with Diagonal/I
    //// Make sure ORDER is correct
    ////// (determine if want traits w/in individs or individs w/in traits)
    KRinv = cs_kroneckerI(Rinv, dimZWG[2]); 
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
    tyRinvy = 0.0;   
    for(k = 0; k < ny[0]; k++){
      tyRinvy += tyRinv->x[k] * y[k];  
    }
    



    // Components of Meyer 1989 eqn 2
    tWKRinv = cs_multiply(tW, KRinv);
    tWKRinvW = cs_multiply(tWKRinv, W);
    // Next creates RHS
    //// Meyer '97 eqn 11
    //// (Note different order of RHS from Meyer '89 eqn 6; Meyer '91 eqn 4)
    cs_gaxpy(tWKRinv, y, RHS->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x` 

  }    // end if lambda=FALSE

  // Now take transpose of transpose to correctly order (don't ask why)
  ttWKRinvW = cs_transpose(tWKRinvW, true);
  cs_spfree(tWKRinvW);
  tWKRinvW = cs_transpose(ttWKRinvW, true);
  cs_spfree(ttWKRinvW);





  // 1c Now make coefficient matrix of MME
  //// form Kronecker products for each G and ginverse element (i.e., I or geninv)
  for(g = 0; g < nG; g++){
    if(ndgeninv[g] == 0){   // Diagonal matrix kronecker product: Ginv %x% I
      KGinv[g] = cs_kroneckerI(Ginv[g], dimZWG[4+g]);
    }else{  // generalized inverse kronecker product: Ginv %x% geninv
      KGinv[g] = cs_kroneckerA(Ginv[g], geninv[g]);
    }
  }
  //// Construct C
  if(nG > 0){
    Ctmp = cs_omega(KGinv, nG, Bpinv);
    C = cs_add(tWKRinvW, Ctmp, 1.0, 1.0);
  }else{
    C = cs_add(tWKRinvW, Bpinv, 1.0, 1.0);
  }




if(v[0] > 3){
  took = toc(t); 
  Rprintf("%6.4f sec.: initial cpp setup (to get C)\n", took);
  t = tic();
}






  //1d Find best order of C pivots!
  //// Graser et al. 1987 (p1363) when singular C, |C| not invariant to order
  //// of C, therefore same order (of pivots) must be used in each loglik iteration
  //// Hadfield (2010, MCMCglmm paper App B): details on chol, ordering, and updating
  
  //// TODO Supernodal decomposition??? Can do it with Csparse? Need something else?
  // Symbolic Cholesky factorization of C
  sLc = cs_schol(1, C);
  if(sLc == NULL){
    error("FAILED: symbolic Cholesky factorization of Coefficient matrix (`C`)\n");
  }




if(v[0] > 3){
  took = toc(t);
  Rprintf("  %6.4f sec.: initial cpp cs_schol(C)\n", took);
  t = tic();
}





  // Create variables needed for log-likelihood calculations
  // Used for log(|R|) and log(|G|) <-- Meyer 1989 (univar.) & 1991 (multivar.)
  for(g = 0; g < nG; g++){
    rfxlvls[g] = dimZgs[2*g+1];
  }
  nminfrfx = nminffx[0] - nr; 









if(v[0] > 3){
  took = toc(t);
  Rprintf("  %6.4f sec.: rest of initial cpp setup\n", took);
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
    if(i == 0){vitout = 0;}else{vitout = (i+1)%vit[0];}  // always do first iteration
    if(v[0] > 0 && vitout == 0){
      Rprintf("  %i of max %i\n", i+1, maxit[0]);//TODO TIME of DAY format as remlIt
    }


    // For ALL iterations i=1 and greater
    if(i > 0){
      cs_spfree(R); cs_spfree(Rinv);
      for(g = 0; g < nG; g++){
        cs_spfree(G[g]);
        cs_spfree(Ginv[g]);
      }
      cs_spfree(C);
      cs_nfree(Lc);



      // setup R matrix
      //// Assumes, just 1 R matrix 
      si = 0; for(g = 0; g < nG; g++) si += nnzGRs[g];
      dimM = dimGRs[nG];
      R = cs_spalloc(dimM, dimM, nnzGRs[nG], true, false);
      for(k = 0; k < nnzGRs[nG]; k++){
        R->i[k] = iGRs[si+k];
        R->x[k] = nu[si+k];
      }
      for(k = 0; k <= dimM; k++){
        R->p[k] = k*dimM;
      }

      Rinv = cs_inv(R);
/* 
could do a check to make sure R was inverted correctly:
  if(Rinv == NULL){
    error("R-structure %i starting value(s) is/are improperly specified: check that all eigenvalues (`eigen(R)$values`) > 0 and that the cholesky decomposition can be formed (`chol(R)`)\n", nR);
  }
*/



      // setup G matrices
      si = 0;
      for(g = 0; g < nG; g++){
        dimM = dimGRs[g];
        G[g] = cs_spalloc(dimM, dimM, nnzGRs[g], true, false);
        for (k = 0; k < nnzGRs[g]; k++){        
          G[g]->i[k] = iGRs[si+k];
          G[g]->x[k] = nu[si+k];
        }
        si += nnzGRs[g];
        for(k = 0; k<= dimM; k++){
          G[g]->p[k] = k*dimM;
        }
        Ginv[g] = cs_inv(G[g]);
      }


      //////////////////////////////////////////////////////////////////////////
      // 1 Setup to create coeeficient matrix of MME (C)
      //// only do this part for NOT lambda: straight (co)variances
      ////// same terms under lambda model don't contain Rinv so don't change
      if(lambda[0] == 0){
        cs_spfree(tWKRinv);
        cs_spfree(tWKRinvW);

        // Update `KRinv` (Rinv Kronecker with Diagonal/I)
        //// Make sure ORDER is correct
        ////// (determine if want traits w/in individs or individs w/in traits)
        cs_kroneckerIupdate(Rinv, dimZWG[2], KRinv);
        // Reset tyRinv to 0.0 so cs_gaxpy works right
        for(k = 0; k < dimZWG[2]; k++) tyRinv->x[k] = 0.0;
        // ASSUME cs_gaxpy does KRinv %*% y, giving correct values for tyRinv@x
        cs_gaxpy(KRinv, y, tyRinv->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x`
 
        tyRinvy = 0.0;   
        for(k = 0; k < ny[0]; k++){
          tyRinvy += tyRinv->x[k] * y[k];  
        }

        // Components of Meyer 1989 eqn 2
        tWKRinv = cs_multiply(tW, KRinv);
        tWKRinvW = cs_multiply(tWKRinv, W);
        // Next creates RHS
        //// Meyer '97 eqn 11
        //// (Note different order of RHS from Meyer '89 eqn 6; Meyer '91 eqn 4)
        // Reset RHS to 0.0 so cs_gaxpy works right
        for(k = 0; k < dimZWG[3]; k++) RHS->x[k] = 0.0;
        cs_gaxpy(tWKRinv, y, RHS->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x` 

      }    // end if lambda=FALSE

      // Now take transpose of transpose to correctly order (don't ask why)
      ttWKRinvW = cs_transpose(tWKRinvW, true);
      cs_spfree(tWKRinvW);
      tWKRinvW = cs_transpose(ttWKRinvW, true);
      cs_spfree(ttWKRinvW);




      // 1c Now make coefficient matrix of MME
      //// form Kronecker products for each G and ginverse element (i.e., I or geninv)
      for(g = 0; g < nG; g++){
        if(ndgeninv[g] == 0){   // Diagonal matrix kronecker product: Ginv %x% I
          cs_kroneckerIupdate(Ginv[g], dimZWG[4+g], KGinv[g]);
        }else{  // generalized inverse kronecker product: Ginv %x% geninv
          cs_kroneckerAupdate(Ginv[g], geninv[g], KGinv[g]);
        }
      }
      //// Construct C
      if(nG > 0){
        cs_omegaupdate(KGinv, nG, Bpinv, Ctmp); 
        C = cs_add(tWKRinvW, Ctmp, 1.0, 1.0);
      }else{
        C = cs_add(tWKRinvW, Bpinv, 1.0, 1.0);
      }


    }  // End if i>0





if(v[0] > 3){
  took = toc(t);
  Rprintf("  %6.4f sec.: cpp REML i=%i setup\n", took, i);
  t = tic();
}






    // Update Cholesky factorization of C
    Lc = cs_chol(C, sLc);  //TODO update when i>0? (cs_updown)
    if(Lc == NULL){
      error("Mixed Model Coefficient matrix (C) singular: possibly caused by a bad combination of G and R (co)variance parameters\n");
    }





if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec.: cpp REML i=%i cs_chol(C)\n", took, i);
  t = tic();
}








    ////////////////////////////////////////////////////////////////////////////
    // solve MME for BLUEs/BLUPs
    //// Do now, because need solutions as part of log-likelihood calc. (for tyPy)
    //XXX tmpBLUXs must have EXPLICIT 0s, else triangular solve = incorrect answer
    cs_ipvec(sLc->pinv, RHS->x, tmpBLUXs->x, C->n) ;   /* x = P*b */
    cs_lsolve(Lc->L, tmpBLUXs->x) ;                    /* x = L\x */
    cs_ltsolve(Lc->L, tmpBLUXs->x) ;                   /* x = L'\x */
    cs_pvec(sLc->pinv, tmpBLUXs->x, BLUXs->x, C->n) ;  /* b = P'*x */


if(v[0] > 3){
  took = toc(t);
  Rprintf("\t    %6.4f sec.: cpp REML i=%i sln forward/back solve with chol(C)\n", took, i);
  t = tic();
}



    // calculate residuals as r = y - W %*% sln
    // put y in r to be over-written
    for(k = 0; k < ny[0]; k++){
      res[k] = -1.0 * y[k];
    }
    //// below rearranges to -r = W %*% sln - y
    cs_gaxpy(W, BLUXs->x, res);   //  y = A*x+y
    for(k = 0; k < ny[0]; k++){
      res[k] *= -1.0;
    }




if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec.: cpp REML i=%i sln/r calc.\n", took, i);
  t = tic();
}









    ////////////////////////////////////////////////////////////////////////////
    // 5 record log-likelihood, check convergence, & determine next varcomps  
    //// 5a determine log(|C|) and y'Py
    // Boldman and Van Vleck 1991 eqn. 6 (tyPy) and 7 (log|C|)    
    // Meyer 1997 eqn. 13 (tyPy)

    // set diagonal sum of Lc, tyPy, and loglik back to 0/starting values
    tyPy = tyRinvy;
    dsLc = 0.0; loglik = 0.0;  
    // 'by hand' calculate `tyPy = tyRinvy -  crossprod(BLUXs, RHS)`
    for(k = 0; k < dimZWG[3]; k++){
      tyPy -= BLUXs->x[k] * RHS->x[k];
    }
    
    // calculate logDetC from diagonals of Lc
    for(k = 0; k < C->n; k++){
      rw = Lc->L->p[k];
      dsLc += log(Lc->L->x[rw]);
    }
    logDetC = 2.0 * dsLc;

    // Factored out residual variance (only for lambda scale)
    if(lambda[0] == 1) sigma2e = tyPy / nminffx[0]; else sigma2e = 1.0;

    // Construct the log-likelihood (Meyer 1997 eqn. 8)
    //// (firt put together as `-2 log-likelihood`)
    // `log(|R|)`
    if(lambda[0] == 1){
      loglik += nminfrfx * log(sigma2e);
    }else{
      loglik += ny[0] * log(nu[p[0]-1]);  //FIXME won't work when R (co)variances
    }
    if(v[0] > 3) Rprintf("\t log|R|=%6.4f\n", loglik);

    // `log(|G|)`
    //// Meyer 1997 eqn. 9 for lambda=FALSE equation
    //// NOTE: setting `sigma2e=1.0` above when lambda=FALSE means 1 line below
    //// FIXME only works for independent random effects
    if(nG > 0){
      for(g = 0; g < nG; g++){
        //FIXME below assumes only one entry in a G[g]
        loglik += rfxlvls[g] * log(G[g]->x[0] * sigma2e);
      }
    }
      if(v[0] > 3) Rprintf("\t log|G| added=%6.4f\n", loglik);
    loglik += rfxlL[0];    // rfxIncContrib2loglik
      if(v[0] > 3) Rprintf("\t\t log|G| rfxlL(%6.4f) added=%6.4f\n", rfxlL[0], loglik);

    // log(|C|) + tyPy
    //// if lambda=TRUE then nminffx=tyPy/sigma2e simplified below
    ////// (because sigma2e=tyPy/nminffx)
    if(lambda[0] == 1){
      loglik += logDetC + nminffx[0];
    }else{
      loglik += logDetC + tyPy;
    }
      if(v[0] > 3) Rprintf("\t log|C|+tyPy added=%6.4f\n", loglik);
    
    // Multiply by -0.5 to calculate `loglik` from `-2loglik`
    loglik *= -0.5;
      if(v[0] > 3) Rprintf("\t multiplied by -0.5=%6.4f\n", loglik);






if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec.: cpp REML i=%i log-likelihood calc.\n", took, i);
  t = tic();
}




    //XXX	**** END LOG-LIKELIHOOD CALCULATION **** 		XXX
    //XXX	**** END LOG-LIKELIHOOD CALCULATION **** 		XXX
    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/





    /* `itc` elements for each iteration:
    nus | sigma2e | tyPy | logDetC | loglik | time
    */
    for(k = 0; k < p[0]; k++){
      itMat[itc] += nu[k];
      itc++;
    }
    itMat[itc] += sigma2e; itc++;
    itMat[itc] += tyPy;    itc++;
    itMat[itc] += logDetC; itc++;
    itMat[itc] += loglik;  itc++;





if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec.: cpp REML i=%i itMat recording took\n", took, i);
  t = tic();
}



    ////////////////////////////////////////////////////////////////////////////
    // 5c check convergence criteria
    //// Knight ('08 ch. 6): Searle et al.'92 & Longford '93 discuss diff. crit. types
    //// Appendix 2 of WOMBAT help manual for 4 criteria specified
    for(k = 0; k <= 4; k++) cc[k] = 0;     // Keep track of sum as last position
    d = 0.0; cc2 = 0.0; cc2d = 0.0;
    if(i > 0){
      // Change in log-likelihood: Wombat 1
      cc[0] += ((itMat[itc-6-p[0]] - loglik) < cctol[0]);
        cc[4] += cc[0];
      // Change in nu: wombat eqn. A.1 also Knight 2008 eqn. 6.1
      for(k = 0; k < p[0]; k++){
        d = nu[k] - itMat[itc-(9+2*p[0])+k];
        cc2 += d*d;
        cc2d += nu[k] * nu[k];
      }
      cc[1] += (sqrt(cc2 / cc2d) < cctol[1]);
        cc[4] += cc[1];
      // 3 & 4 are only for optimzation algorithms which produce derivatives
      if(algit[i] > 0){    // 0=EM, 1=AI
        // Norm of gradient vector: wombt eqn. A.2
//TODO Does this go here or maybe after AI
//// Does this step happen for last AI matrix (i-1) or current (i)?
        d = 0.0;
        for(k = 0; k < p[0]; k++) d += dLdnu[k] * dLdnu[k];
        cc[2] += (sqrt(d) < cctol[2]);
          cc[4] += cc[2];
        // Newton decrement: wombat eqn A.3 and Boyd & Vandenberghe 2004
          //TODO
          //  cc[4] += cc[3];
      }  // end if AI
    }  // end convergence criteria if i>0




if(v[0] > 3){
  took = toc(t);
  Rprintf("    %6.4f sec.: cpp REML i=%i convergence crit. calc.\n", took, i);
  t = tic();
}



    ////////////////////////////////////////////////////////////////////////////
    // 5d Determine next (co)variance parameters to evaluate: REML NOT CONVERGED
    if(cc[4] < 3){

      /////////////////////////////
      // Expectation Maximization
      /////////////////////////////
      if(algit[i] == 0){
        if(v[0] > 1 && vitout == 0) Rprintf("\tEM to find next theta");
        if(!cs_em(BLUXs, nu, Cinv_ii,
	    nG, rfxlvls, nffx, ndgeninv, geninv, Lc->L, sLc->pinv)){
          error("Unsuccessful EM algorithm in iteration %i\n", i);
        }
        // Calculate EM for residual:
        //// crossprod(y, r) / nminffx
        d = 0.0;
        for(k = 0; k < ny[0]; k++) d += y[k] * res[k]; 
        nu[nG] = d / nminffx[0];
      }  // end EM
      /////////////////////////////




      /////////////////////////////
      // Average Information
      /////////////////////////////
      if(algit[i] == 1){
        if(v[0] > 1 && vitout == 0) Rprintf("\tAI to find next theta");
        if(aiformed == 1) cs_spfree(AI);
        if(lambda[0] == 1){
          AI = cs_ai(BLUXs, Ginv, R, 0, 0,
	      y, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	      0, sigma2e);
          if(AI == NULL) error("Unsuccessful AI algorithm in iteration %i\n", i);

if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate AI", took);
  t = tic();
}
       
          if(!cs_gradFun(nu, dLdnu, Cinv_ii,
	      ny[0], p[0], nG, rfxlvls, nffx, ndgeninv,
	      geninv, BLUXs, Lc->L, sLc->pinv, 
              sigma2e,    // 1.0 if lambda=FALSE
	      0, res)){      // 0 if lambda=TRUE
	      
            error("Unsuccessful gradient calculation in iteration %i\n", i);

if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  t = tic();
}
          }  // end if cs_gradFun

        }else{
          AI = cs_ai(BLUXs, Ginv, R, KRinv, tWKRinv,
	      res, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	      nG, 1.0);
          if(AI == NULL) error("Unsuccessful AI algorithm in iteration %i\n", i);

 if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate AI", took);
  t = tic();
}
          
          if(!cs_gradFun(nu, dLdnu, Cinv_ii,
	      ny[0], p[0], nG, rfxlvls, nffx, ndgeninv,
	      geninv, BLUXs, Lc->L, sLc->pinv, 
              1.0,    // 1.0 if lambda=FALSE
	      nG, res)){      // 0 if lambda=TRUE

            error("Unsuccessful gradient calculation in iteration %i\n", i);

if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  t = tic();
}

          }  // end if cs_gradFun
        }  // end if/else lambda

        //TODO do I need to check convergence criteria here (i.e., cc[3:4])


        // Find next set of parameters using a quasi-Newton method/algorithm
        //// Meyer 1989 pp. 326-327 describes quasi-Newton methods 
//TODO see Meyer 1997 eqn 58 for Marquardt 1963: theta_t+1=theta_t - (H_t + k_t * I)^{-1} g_t 
        //// Mrode 2005 eqn 11.4
        //// Johnson and Thompson 1995 eqn 12
        //////(though gremlin uses `+` instead of J & T '95 `-` because
        ////// gremlin multiplies gradient by -0.5 in `gradFun()`)

        //// Check/modify AI matrix to 'ensure' positive definiteness
        //// `fI` is factor to adjust AI matrix
        ////// (e.g., Meyer 1997 eqn 58 and WOMBAT manual A.5 strategy 3b)

/*
        check if any small/zero eigenvalues of AI
        AIeigvals <- eigen(AI)$values
          d <- (3*10^-6) * AIeigvals[1]
          f <- max(0, d - AIeigvals[nrow(AI)])
        fI <- f * diag(x = 1, nrow = nrow(AI))
	//////// modified 'Hessian'
        H <- fI + AI
*/
        // Check if AI can be inverted
        if(aiformed == 1){
          cs_sfree(sLai);  // each time in case AI pattern changes (fix components)
          cs_nfree(Lai);
        }
        sLai = cs_schol(1, AI);
        Lai = cs_chol(AI, sLai);
        aiformed = 1;
        if(Lai == NULL){
          if(v[0] > 1){
            Rprintf("AI cholesky decomposition failed:\n\t AI matrix may be singular - switching to an iteration of the EM algorithm\n");
          } // end if v>1

          //// if AI cannot be inverted do EM
          //////////////  TEMPORARY EM    /////////
          if(v[0] > 1 && vitout == 0) Rprintf("\t\tEM to find next theta");
          if(!cs_em(BLUXs, nu, Cinv_ii,
	      nG, rfxlvls, nffx, ndgeninv, geninv, Lc->L, sLc->pinv)){
            error("Unsuccessful EM algorithm in iteration %i\n", i);
          }
          // Calculate EM for residual:
          //// crossprod(y, r) / nminffx
          d = 0.0;
          for(k = 0; k < ny[0]; k++) d += y[k] * res[k]; 
          nu[nG] = d / nminffx[0];
          algit[i] = 0;  // switch algorithm to EM so itMat of output is accurate
        // end EM
        } else{  //<-- end if AI cannot be inverted


          AIinv = cs_inv(AI);
          //Hinv = cs_inv(H);
//TODO need a check that not proposing negative/0 variance or |correlation|>1
//// Require restraining naughty components
          // fill `nu` with parameters proposed for next iteration
          ////  cs_gaxpy is y = A*x+y where y=nu and x=dLdnu
          cs_gaxpy(AIinv, dLdnu, nu);

          for(g = 0; g < p[0]; g++){
            //FIXME check variances and cov/corr separately
            if(nu[g] < ezero[0]){
              if(v[0] > 1){
                Rprintf("Variance component %i fixed to zero\n", g+1);
              }
              nu[g] = ezero[0];  //FIXME TODO!!!
            }  // end if nu < ezero
          }  // end for g

        }  //<-- end else AI can be inverted
        cs_spfree(AIinv);
      }  // end AI
      /////////////////////////////











if(v[0] > 3){
  took = toc(t);
  Rprintf(": %6.4f sec.\n", took, i);
  t = tic();
}



    }  // end if REML did not converge












    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if(lambda[0] == 1) nu[nG] = 1.0;  // keep R=1 (R factored out)
                                      // TODO FIXME if >1 residual (co)variance

    took = toc(T);                 // Capture cpu clock time for i REML iteration
    // V=1 LEVEL of OUTPUT
    if(v[0] > 0 && vitout == 0){ 
      Rprintf("\n\tlL:%6.6f", loglik);
      Rprintf("\t\ttook %6.4f sec.\n", took); //TODO format units if >60 (and do for all Rprintf(took))
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
          if(algit[i] == 1){
            Rprintf("\tgradient | AI\n");
            Rprintf("\t-------- |--------\n");
            for(g = 0; g < p[0]; g++){
              Rprintf("\t%6.4g\t| ", dLdnu[g]);  // prints gradient value
              // print AI[g, ]
              for(k = 0; k < p[0]; k++){
                for(si = AI->p[k]; si < AI->p[k+1]; si++){
                  // prints AI[g, k]
                  if(AI->i[si] == g) Rprintf(" %6.4f", AI->x[si]);
                }  // end for si (rows of AI column)
              }  // end for k (column of AI)
              Rprintf("\n");  // start next row of output/AI
            }  // end for g (gth parameter/row of AI)  
          }  // end algit == AI
        }  // end if v>2
      }  // end if v>1
    }  // end if v>0


    itMat[itc] += round(took*10) / 10;              // gives 1 decimal place
    itc++;
 



   
    // Determine if model has converged
    //FIXME: change number of parameters that must be true as add criteria
    if(cc[4] > 2){ // FIXME: change number of parameters that must be true
      Rprintf("\n\nREML converged\n\n");      
      break;
    }

  }  // end i for loop 
  maxit[0] = i;    // record number of iterations

  //////////////////////////////////////////////////////////////////////////////
  //XXX									     XXX
  /*		END reml iterations		END reml iterations	      */
  //XXX									     XXX
  //////////////////////////////////////////////////////////////////////////////



  // Calculate Cinv_ii and AI
  cs_chol2inv_ii(Lc->L, sLc->pinv, Cinv_ii, 0);
  //// Average Information
  if(aiformed == 1) cs_spfree(AI);
  if(lambda[0] == 1){
    AI = cs_ai(BLUXs, Ginv, R, 0, 0,
	  y, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	  0, sigma2e);
    if(AI == NULL) error("Unsuccessful AI algorithm at convergence %i\n", i);

if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate AI", took);
  t = tic();
}
       
/*
    if(!cs_gradFun(nu, dLdnu, Cinv_ii,
   	    ny[0], p[0], nG, rfxlvls, nffx, ndgeninv,
	    geninv, BLUXs, Lc->L, sLc->pinv, 
            sigma2e,    // 1.0 if lambda=FALSE
	    0, res)){      // 0 if lambda=TRUE
	      
    error("Unsuccessful gradient calculation  at convergence %i\n", i);

if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  t = tic();
}
    }  // end if cs_gradFun
*/

  }else{
    AI = cs_ai(BLUXs, Ginv, R, KRinv, tWKRinv,
          res, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	  nG, 1.0);
    if(AI == NULL) error("Unsuccessful AI algorithm  at convergence %i\n", i);

 if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate AI", took);
  t = tic();
}

/*          
    if(!cs_gradFun(nu, dLdnu, Cinv_ii,
          ny[0], p[0], nG, rfxlvls, nffx, ndgeninv,
	  geninv, BLUXs, Lc->L, sLc->pinv, 
          1.0,    // 1.0 if lambda=FALSE
	  nG, res)){      // 0 if lambda=TRUE

    error("Unsuccessful gradient calculation  at convergence %i\n", i);

if(v[0] > 3){
  took = toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  t = tic();
}

    }  // end if cs_gradFun
*/
  }  // end if/else lambda



if(v[0] > 3) t = tic();

  // Pass information to R
  //// Solution vector
  for(k = 0; k < dimZWG[3]; k++) sln[k] += BLUXs->x[k];


//XXX keep this commented out bit for ref on how to un-pack a cs matrix
/* No longer needed - cs_em automatically replaces and don't make Cinv
  //// diagonals of Cinv (sln sampling variances) 
  for(k = 0; k < Cn; k++){
    for(i = Cinv->p[k]; i < Cinv->p[k+1]; i++){
      if(Cinv->i[i] == k){
        Cinv_ii[k] += Cinv->x[i];
        break; 
      }  // end if
    }  // end for i
  }  // end for k
*/

  //// return AI
  if(CS_CSC(AI)){
    for(k = 0; k < p[0]; k++){
      for(g = AI->p[k]; g < AI->p[k+1]; g++){
        i = AI->i[g];
        AIvec[k*p[0]+i] += AI->x[g];
      }
    }
    cs_spfree(AI);
  }  // end if AI NOT NULL


  // return permutation matrix of symbolic Cholesky factorization of C
  for(k = 0; k < C->m; k++) sLcPinv[k] += sLc->pinv[k];


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////











  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  cs_spfree(Bpinv);
  cs_spfree(W); cs_spfree(tW);
  cs_spfree(R); cs_spfree(Rinv); cs_spfree(tyRinv);
  if(lambda[0] == 0){
    cs_spfree(KRinv); cs_spfree(tWKRinv); 
  }
  cs_spfree(tWKRinvW);
  cs_spfree(Ctmp); cs_spfree(C);
  cs_spfree(RHS); cs_spfree(tmpBLUXs); cs_spfree(BLUXs);

  cs_sfree(sLc);
  cs_nfree(Lc);

  if(aiformed == 1){
    cs_sfree(sLai);
    cs_nfree(Lai);
  }

//
  for(g = 0; g < nG; g++){
    if(ndgeninv[g] == 1) cs_spfree(geninv[g]);
    cs_spfree(G[g]); cs_spfree(Ginv[g]);
    cs_spfree(KGinv[g]);
  } 
//
  delete [] geninv;
  delete [] G; delete [] Ginv;
  delete [] KGinv;
//
  delete [] rfxlvls;
  delete [] cc;


if(v[0] > 3){
  took = toc(t);
  Rprintf("%6.4f sec.: cpp post-REML freeing-up\n", took);
}

}
}

