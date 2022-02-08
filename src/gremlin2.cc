#include "gremlincc.h"



extern "C"{  
void ugremlin2(
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
/*23*/  int *conv,		// constraint codes (F=0)
/*24*/  double *bound,          // boundaries (0:p-1=LB | p:2p-1=UB) Fixed/NA=0
/*25*/	int *nnzBpinv,		// inverse Fixed effects prior matrix (no prior=0s)
/*26*/	int *iBpinv,
/*27*/	int *pBpinv,
/*28*/	double *xBpinv,
/*29*/	double *dLdnu,		// gradient vector (1st deriv. of lL / nus)
/*30*/	double *AIvec,		// column-wise vector of average information matrix (2nd deriv. of lL / nu)
/*31*/	double *sln,		// solution vector
/*32*/  double *Cinv_ii,	// diagonal of Cinv: sampling variances of BLUXs
/*33*/	double *res,		// residual vector
/*34*/	double *itMat,		// parameter information at each iteration
/*35*/	int *algit,		// algorithm for optimization at each iteration
/*36*/  int *fdit,		// Finite Difference algorithm at each iteration
/*37*/	int *maxit,		// maximum No. iterations
/*38*/  double *step,		// initial/default step-halving value
/*39*/	double *cctol,		// convergence criteria tolerances
/*40*/	double *ezero,		// effective zero
/*41*/  double *einf,           // effective +/- max values
/*42*/	int *v,			// verbosity level of  output messages
/*43*/	int *vit,		// at what iterations to output messages
/*44*/  int *sLcPinv		// empty Cholesky permutation vector
		
){

  cs	 *Bpinv, *W, *tW,
	 *R, *Rinv, *KRinv, *tWKRinv, *tWKRinvW, *ttWKRinvW,
	 *Ctmp, *C,
	 *RHS, *tmpBLUXs, *BLUXs,
	 *AI;
	 //TODO when implement Newton Decrement method need *invAI;

  css    *sLc;
  csn    *Lc;

  int    nG = nGR[0];

  cs*    *geninv = new cs*[nG];
  cs*    *G = new cs*[nG];
  cs*    *Ginv = new cs*[nG];
  cs* 	 *KGinv = new cs*[nG];

  double t[2], T[2], took, tyPy, tyRinvy, logDetC, sigma2e,
         loglik, d, cc2, cc2d, f, stpVal;

  int 	 g, i, k, si, si2, vitout,
	 itc = 0,
         dimM,      // GENERIC dimension of a matrix variable to be REUSED 
	 nffx, nminfrfx,
	 nr = dimZWG[1],
	 bd;  

  int	 *rfxlvls = new int[nG];

  int	 *cc = new int[5];

  int    *con = new int[p[0]];
    for(k = 0; k < p[0]; k++) con[k] = conv[k];  // to retain original values

  // first p[0]=Lwr. Bounds, second p[0]=Upr. Bounds
  int    *wchBd = new int[2*p[0]];  
    for(k = 0; k < 2*p[0]; k++) wchBd[k] = 0;  // initialize entire vector

  double  *w = new double[dimZWG[3]];

  g = (lambda[0]) ? nG : nG+1;
    double  *tugug = new double[g];  // includes crossprod(residual) when !lambda

  double  *trace = new double[nG];

  double  *newnu = new double[p[0]];

  double  *dnu = new double[p[0]];

if(v[0] > 3) simple_tic(t);







//FIXME  FIXME    FIXME   DELETEME    DELETEME  FIXME
  Rprintf("\n\n\t\t**GREMLIN 2** called\n\n");
//XXX	XXX	XXX	XXX XXX			XXX




  sigma2e = (lambda[0]) ? 0.0 : 1.0;

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
    tmpBLUXs = cs_spalloc(dimZWG[3], 1, dimZWG[3], true, false);  // sln matrix
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


  // Create variables needed for log-likelihood calculations
  // Used for log(|R|) and log(|G|) <-- Meyer 1989 (univar.) & 1991 (multivar.)
  for(g = 0; g < nG; g++){
    rfxlvls[g] = dimZgs[2*g+1];
  }
  // only use `nminfrfx` when lambda=TRUE
  //// but just do simple calculation regardless, to avoid if statement
  nminfrfx = nminffx[0] - nr; 

  // initialize to zero: will change for lambda=TRUE else always be passed to reml
  //// as 0.0 to be calculated within cs_reml() for each G/R
  tyRinvy = 0.0;
  if(lambda[0] == 1){
    // below don't change between REML iterations for lambda=TRUE
    tWKRinvW = cs_multiply(tW, W);
    // Fill in RHS
    cs_gaxpy(tW, y, RHS->x);  // y = A*x+y XXX my `y` is cs_gaxpy's `x` 
    // form tyy (crossprod(y)), but call it tyRinvy for conciseness
    //// either non-zero for all iterations or 0.0 passed in when lambda=FALSE
    for(k = 0; k < ny[0]; k++){
      tyRinvy += y[k] * y[k];
    }
  }  // end if lambda=TRUE


  // setup empty AI matrix
  AI = cs_spalloc(p[0], p[0], p[0]*p[0], true, false);



if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("  %6.4f sec.: initial cpp setup\n", took);
}








  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //XXX									   XXX
  //XXX									   XXX
  //XXX									   XXX
  /*		REML ITERATIONS			REML ITERATIONS		    */
  //XXX									   XXX
  //XXX									   XXX
  //XXX									   XXX
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  for(i = 0; i < maxit[0]; i++){
    simple_tic(T);
    if(i == 0){vitout = 0;}else{vitout = (i+1)%vit[0];}  // always do first iteration
    if(v[0] > 0 && vitout == 0){
      Rprintf("  %i of max %i\t", i+1, maxit[0]);//TODO TIME of DAY format as remlIt
      // V=2 LEVEL of OUTPUT
      if(v[0] > 1){
        Rprintf("\n");
        // output names of variables
        for(g = 0; g < p[0]; g++) Rprintf("\tV%i", g+1);
        Rprintf("\n");
        // output underlines
        //// =No. characters (2 for V1-V9 and 3 for V10 and up)
        for(g = 0; g < p[0]; g++){
          if(g < 10) Rprintf("\t--"); else Rprintf("\t---");  
        }
        Rprintf("\n");
        // output variable values
        for(g = 0; g < p[0]; g++) Rprintf("\t%6.4g", nu[g]);
        Rprintf("\n");
      }  // end v>1
    }  // end v>0


    // For ALL iterations i=1 and greater
    if(i > 0){
      cs_spfree(R); cs_spfree(Rinv);
      for(g = 0; g < nG; g++){
        cs_spfree(G[g]);
        cs_spfree(Ginv[g]);
      }
      cs_nfree(Lc);
      
      if(lambda[0] == 0){
        cs_spfree(tWKRinv);
        cs_spfree(tWKRinvW);
      }
      
    }  // End if i>0

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



    // need to set up sLc once and not again
    //// just construct matrices needed to make C
    if(lambda[0] == 0){
      // NOT lambda: Straight (co)variances
      // Rinv Kronecker with Diagonal/I
      //// Make sure ORDER is correct
      ////// (determine if want traits w/in individs or individs w/in traits)
      if(i == 0){
        KRinv = cs_kroneckerI(Rinv, dimZWG[2]); 
      } else cs_kroneckerIupdate(Rinv, dimZWG[2], KRinv);

      // Components of Meyer 1989 eqn 2
      tWKRinv = cs_multiply(tW, KRinv);
      tWKRinvW = cs_multiply(tWKRinv, W);
    }
     
    // Now take transpose of transpose to correctly order (don't ask why)
    ttWKRinvW = cs_transpose(tWKRinvW, true);
    cs_spfree(tWKRinvW);
    tWKRinvW = cs_transpose(ttWKRinvW, true);
    cs_spfree(ttWKRinvW);
    
    
    // 1c Now make coefficient matrix of MME
    if(i == 0){   
      // form Kronecker products for each G and ginverse element (i.e., I or geninv)
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
  took = simple_toc(t); 
  Rprintf("%6.4f sec.: initial cpp setup (to get C)\n", took);
  simple_tic(t);
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
      cs_spfree(C);  // make inside cs_reml so don't need again here 

if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("  %6.4f sec.: initial cpp cs_schol(C)\n", took);
  simple_tic(t); 
}


    
    } // end setup if first iteration (i == 0)  
    


    Lc = cs_reml(ny[0], dimZWG, nG, p[0], y,
        Bpinv, W, tW, rfxlvls, rfxlL[0],
        R, Rinv, G, Ginv, ndgeninv, geninv,
        KRinv, KGinv, tWKRinv, tWKRinvW, Ctmp,
        RHS, tmpBLUXs, BLUXs, res,
        sLc, 
        &tyPy, &logDetC, &sigma2e,
        tyRinvy,
        nminffx[0],
        &loglik,
        i, v[0], vitout, lambda[0]); 
    if(loglik == 0.0){
      error("\nUnsuccessful REML calculation: iteration %i", i);
    }




    //XXX	**** END LOG-LIKELIHOOD CALCULATION **** 		XXX
    /**************************************************************************/


    simple_tic(t);
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
  took = simple_toc(t);
  Rprintf("\n    %6.4f sec.: cpp REML i=%i itMat recording", took, i);
  simple_tic(t);
}

    // V=1 LEVEL of OUTPUT
    //// same line if minimal output, different lines if more extensive output
    if(v[0] > 1 && vitout == 0) Rprintf("\n");  
    if(v[0] > 0 && vitout == 0) Rprintf("\tlL:%6.6f", loglik);



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
    }  // end convergence criteria if i>0







    ////////////////////////////////////////////////////////////////////////////
    // 5d Determine next (co)variance parameters to evaluate:

    /////////////////////////////
    // Average Information
    /////////////////////////////
    if(algit[i] > 0){
      if(v[0] > 1 && vitout == 0) Rprintf("\n\tAI to find next nu");
      f = 0.0;  // initialize to a default/no alteration of Hessian


      // Gradient via Analytical (vs. Finite Difference)
      if(algit[i] == 1){  

        if(!tugugFun(tugug, w, nG, rfxlvls, con,
	    nffx, ndgeninv, geninv, BLUXs)){
          error("\nUnsuccessful tugug calculation: AI algorithm in iteration %i", i);
        }
if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.6f sec.: calculate tugug(s)", took);
  simple_tic(t);
}

        if(!traceFun(trace, w, nG, rfxlvls,
	    nffx, ndgeninv, geninv, BLUXs->m, Lc->L, sLc->pinv)){
          error("\nUnsuccessful trace calculation: AI algorithm in iteration %i", i);
        }
if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.6f sec.: calculate trace(s)", took);
  simple_tic(t);
}
      }  // end if gradient via analytical (instead of finite difference)
      
      
           
      if(lambda[0] == 1){
        cs_spfree(AI);  //TODO how pass if not initialized
        AI = cs_ai(BLUXs, Ginv, R, 0, 0,
	      y, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	      0, sigma2e);
	if(AI == 0){
          error("\nUnsuccessful AI algorithm iteration %i", i);
        }
if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.6f sec.: calculate AI matrix", took);
  simple_tic(t);
}


        if(algit[i] == 1){
          if(!cs_gradFun(nu, dLdnu,
	      tugug, trace, con,
	      ny[0], nG, rfxlvls, nffx,
              sigma2e,    // 1.0 if lambda=FALSE
	      0, res)){      // 0 if lambda=TRUE
	      
            error("\nUnsuccessful gradient calculation iteration %i", i);
          }  // end if cs_gradFun
if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  simple_tic(t);
}
        }  // end if algit=1
               

      } else{
          // when lambda = FALSE
          cs_spfree(AI);
          AI = cs_ai(BLUXs, Ginv, R, KRinv, tWKRinv,
	      res, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	      nG, 1.0);
	  if(AI == NULL){
          error("\nUnsuccessful AI algorithm in iteration %i", i);
          }

if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.6f sec.: calculate AI matrix", took);
  simple_tic(t);
}


          // algit 1 is AI with analytical gradient calculation
          if(algit[i] == 1){
            if(!cs_gradFun(nu, dLdnu,
	      tugug, trace, con,
	      ny[0], nG, rfxlvls, nffx,
              1.0,    // 1.0 if lambda=FALSE
	      nG, res)){      // 0 if lambda=TRUE
            error("\nUnsuccessful gradient calculation iteration %i", i);
          }  // end if cs_gradFun

if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  simple_tic(t);
}
          }  // end if algit=1     

        }  // end if/else lambda





      // algit 2 is AI with finite differences for the gradient calculation
      //// specify 1 way regardless of lambda
      if(algit[i] == 2){
        if(!cs_gradFun_fd(nu, fdit[i], dLdnu, loglik, con,
            ny[0], dimZWG, nG, p[0], y,
            Bpinv, W, tW, rfxlvls, rfxlL[0],
            ndgeninv, geninv, KRinv,
            Ctmp, RHS, tmpBLUXs, BLUXs,
            sLc,
            tyRinvy,
            nminffx[0],
            nnzGRs, dimGRs, iGRs, lambda[0])){      // 0 if lambda=TRUE
          error("\nUnsuccessful finite difference gradient calculation iteration %i", i);
        }  // end if cs_gradFun

if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate gradient", took);
  simple_tic(t);
}
      }  // end if algit=2




      // Find next set of parameters using a quasi-Newton method/algorithm
      //// Meyer 1989 pp. 326-327 describes quasi-Newton methods 
//TODO see Meyer 1997 eqn 58 for Marquardt 1963: theta_t+1=theta_t - (H_t + k_t * I)^{-1} g_t 
      // What I do below is similar: except k_t=f
      //// Mrode 2005 eqn 11.4
      //// Johnson and Thompson 1995 eqn 12
      //////(though gremlin uses `+` instead of J & T '95 `-` because
      ////// gremlin multiplies gradient by -0.5 in `gradFun()`)

      for(k = 0; k < p[0]; k++) newnu[k] = 0.0;  // clear

      if(!qNewtRhap(nu, newnu, dLdnu, AI,
		p[0], con, wchBd, f, ezero, v[0])){
        // AI algorithm failed: do EM
        if(v[0] > 1 && vitout == 0) Rprintf("\n\t\tAI failed, switching to EM");
        algit[i] = 0;  // switch algorithm to EM
      } 


      // if AI working so far
      // calculate `dnu` - proposed change in nu parameters (`Hinv %*% grad`)
      // First, implement step-reduction (if necessary)
      //// Rule: if proposed values >200% change in any parameter
      //// then implement step reduction (`step[0]` default) else do not

      si2 = 0;  // Keep track (avoid infinite loops): 
                //// get as many tries as numbers of parameters
      bd = 1;  // initialize so enters while loop
      while(algit[i] > 0 && bd == 1 && si2 < p[0]){
        si = 0;
        stpVal = 1.0;
        for(k = 0; k < p[0]; k++){
          if(con[k] == 0 || con[k] == 3) continue;
          dnu[si] = newnu[si] - nu[k];
          if(abs(dnu[si] / nu[k]) > 2.0){
            stpVal = step[0];
          }
          si++;
        }  // end for k
 
        if(stpVal == step[0]){  // if TRUE then implement step-reduction
          if(v[0] > 2) Rprintf("\n\tstep reduction: %6.4f", stpVal);
          si = 0;
          for(k = 0; k < p[0]; k++){
            if(con[k] == 0 || con[k] == 3) continue;
            newnu[si] = nu[k] + dnu[si] * stpVal;
            si++;
          }
        }  
    
        // Second, check for indecent proposals
        si = 0; bd = 0;
        for(k = 0; k < 2*p[0]; k++) wchBd[k] = 0;  // reset, esp. b/c was used above
        for(g = 0; g < p[0]; g++){
          if(con[g] == 0) continue;
          if(newnu[si] <= bound[g]){
            bd = 1;
            wchBd[g] += 1;
          }
          if(newnu[si] >= bound[p[0] + g]){
            bd = 1;
            wchBd[p[0] + g] += 1;
          }
          si++;
        }
 
        // restrain naughty components that enacted indecent proposals
        if(bd == 1){
          if(v[0] > 0){
            Rprintf("\n(co)variance component(s) in `thetav` vector (position:");
            for(g = 0; g < p[0]; g++){
              if(wchBd[g] == 1 || wchBd[p[0] + g] == 1) Rprintf(" %i ", g+1);
            }
            Rprintf(") restrained inside boundaries\t");
          }
          si = 0;
          for(g = 0; g < p[0]; g++){
            if(con[g] == 0) continue;
            if(wchBd[g] == 1){
              newnu[si] = bound[g] + ezero[0];
              con[g] = 3;
            }
            if(wchBd[p[0] + g] == 1){
              newnu[si] = bound[p[0] + g] - ezero[0];
              con[g] = 3;
            }
            si++;
          }
    

          // Re-calculate parameter updates, CONDITIONAL on restrained values
          //// Gilmour 2019 AI REML in Practice. J. Anim. Breed. Genet.
          if(!qNewtRhap(nu, newnu, dLdnu, AI,
     		        p[0], con, wchBd, f, ezero, v[0])){ 
            // AI algorithm failed: do EM
            if(v[0] > 1 && vitout == 0) Rprintf("\n\t\tAI failed, switching to EM");
            algit[i] = 0;  // switch algorithm to EM
          } 

 
        }  // end if bd/indecent proposals

        si2++;

      }  // end while AI and bd


      if(algit[i] > 0){
        // CONVERGENCE CRITERIA 3 and 4
        //// Appendix 2 of WOMBAT help manual for 4 criteria specified
        // Norm of gradient vector: wombt eqn. A.2
        d = 0.0;
        for(k = 0; k < p[0]; k++){
          if(con[k] == 0 || con[k] == 3) continue;
          d += dLdnu[k] * dLdnu[k];
        }
        cc[2] += (sqrt(d) < cctol[2]);
          cc[4] += cc[2];
        // Newton decrement: wombat eqn A.3 and Boyd & Vandenberghe 2004
        /* TODO: determine `cctol` value - FIXME truned off for now 
        invAI = cs_inv_withDiagMod(AI, con, wchBd, ezero, v[0]);
        d = 0.0;
        for(g = 0; g < p[0]; g++){
          for(k = invAI->p[g]; k < invAI->p[g+1]; k++){
            d += dLdnu[ invAI->i[k] ] * dLdnu[g] * invAI->x[ invAI->i[k] ];
          }  // end for kth row
        }  // end for gth column 
        cc[3] += (-1 * d) < cctol[3];
          cc[4] += cc[3];
        */

        // unpack newnu into nu
        si = 0;
        for(k = 0; k < p[0]; k++){
          if(con[k] == 0) continue;
          nu[k] = newnu[si];
          si++;
        }


      }  // end if AI did NOT fail


      // Remove any boundary constraint codes from previous iterations
      //// restore to original coded value
      for(g = 0; g < p[0]; g++) if(con[g] == 3) con[g] = conv[g];







    }  // end AI
    ///////////////////////////////
    ///////////////////////////////






    /////////////////////////////
    // Expectation Maximization
    /////////////////////////////
    // Place after AI: if AI fails, change value of algit[i] and will do EM here
    if(algit[i] == 0 && cc[4] < 2){
      if(v[0] > 1 && vitout == 0) Rprintf("\n\tEM to find next nu");

      if(!tugugFun(tugug, w, nG, rfxlvls, con,
	    nffx, ndgeninv, geninv, BLUXs)){
        error("\nUnsuccessful tugug calculation: EM algorithm in iteration %i", i);
      }

      if(!traceFun(trace, w, nG, rfxlvls,
	    nffx, ndgeninv, geninv, BLUXs->m, Lc->L, sLc->pinv)){
        error("\nUnsuccessful trace calculation: EM algorithm in iteration %i", i);
      }

      // calculate EM for G (co)variances:
      //// (tugug + trace ) / qi
      for(g = 0; g < nG; g++){
        if(con[g] == 0) continue;  // skip if parameter is fixed
        nu[g] = (tugug[g] + trace[g] ) / rfxlvls[g];
      }

      // Calculate EM for residual:
      //// crossprod(y, r) / nminffx
      // TODO/FIXME below when >1 Residual
      if(con[nG] != 0){  // only do if Residual is NOT fixed 
        d = 0.0;
        for(k = 0; k < ny[0]; k++) d += y[k] * res[k]; 
        nu[nG] = d / nminffx[0];
      }
    }  // end EM
    /////////////////////////////








    // V=2 LEVEL of OUTPUT
    if(v[0] > 1 && vitout == 0){ 
      // output convergence criteria
      Rprintf("\n\tConvergence crit:");
      for(k = 0; k < 4; k++) Rprintf("%4i", cc[k]);
      Rprintf("\n");
    }  // end v > 1










    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    if(lambda[0] == 1) nu[nG] = 1.0;  // keep R=1 (R factored out)
                                      // TODO FIXME if >1 residual (co)variance

    took = simple_toc(T);        // Capture cpu clock time for i REML iteration
    // V=1 LEVEL of OUTPUT
    if(v[0] > 0 && vitout == 0){
      // V=2 LEVEL of OUTPUT
      if(v[0] > 1 && vitout == 0){ 
        // V=3 LEVEL of OUTPUT
        if(v[0] > 2){
          if(algit[i] > 0){
            // output step-size modification
            Rprintf("\tgradient | AI\n");
            Rprintf("\t-------- |--------\n");
            for(g = 0; g < p[0]; g++){
              Rprintf("\t%6.4g\t| ", dLdnu[g]);  // prints gradient value
              // print AI[g, ]
              for(k = 0; k < p[0]; k++){
                for(si = AI->p[k]; si < AI->p[k+1]; si++){
                  // prints AI[g, k]
                  if(AI->i[si] == g) Rprintf(" %6.4g", AI->x[si]);
                }  // end for si (rows of AI column)
              }  // end for k (column of AI)
              Rprintf("\n");  // start next row of output/AI
            }  // end for g (gth parameter/row of AI)  
          }  // end algit == AI
        }  // end if v > 2
      }  // end if v > 1
      Rprintf("\t\ttook %6.4f sec.\n", took); //TODO format units if >60 (and do for all Rprintf(took))
      // To format units see const *char in: http://www.cplusplus.com/reference/cmath/round/
    }  // end if v > 0


    itMat[itc] += round(took*10) / 10;              // gives 1 decimal place
    itc++;
 



   
    // Determine if model has converged
    //FIXME: change number of parameters that must be true as add criteria 4
    //// (FIXME: also add "bump" give to cc[4] with EM algorithm)
    // if EM algorithm add TRUE for number of extra criteria used in AI
    if(algit[i] == 0) cc[4]++;  
    if(cc[4] > 2){
      if(v[0] > 0) Rprintf("\n***  REML converged  ***\n\n"); 
      i++;     
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
if(v[0] > 3){
  took = simple_toc(t); 
  Rprintf("\n\t    %6.4f sec.: calculate Cinv_ii", took);
  simple_tic(t); 
}

  //// Average Information
  ////// only need to do if did NOT do AI
  if(algit[i] < 1){  
    if(lambda[0] == 1){
      cs_spfree(AI);  //TODO how pass if not initialized
      AI = cs_ai(BLUXs, Ginv, R, 0, 0,
	      y, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	      0, sigma2e);
	if(AI == NULL){
          error("\nUnsuccessful AI algorithm in iteration %i", i);
        }

if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n\t    %6.4f sec.: calculate AI", took);
  simple_tic(t);
}
       

    }else{
      cs_spfree(AI);  //TODO how pass if not initialized
        AI = cs_ai(BLUXs, Ginv, R, KRinv, tWKRinv,
	      res, W, tW, ny[0], p[0], nG, rfxlvls, nffx, Lc->L, sLc->pinv,
	      nG, 1.0);
	if(AI == NULL){
          error("\nUnsuccessful AI algorithm in iteration %i", i);
        }
if(v[0] > 3){
  took = simple_toc(t); 
  Rprintf("\n\t    %6.4f sec.: calculate AI", took);
  simple_tic(t); 
}

    }  // end if/else lambda
  }  // end if algit!=AI


if(v[0] > 3) simple_tic(t);

  // Pass information to R
  //// Solution vector
  for(k = 0; k < dimZWG[3]; k++) sln[k] += BLUXs->x[k];


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

  // return constraint codes from final iteration
  for(k = 0; k < p[0]; k++) conv[k] = con[k];  // to retain original values


  // return permutation matrix of symbolic Cholesky factorization of C
  for(k = 0; k < BLUXs->m; k++) sLcPinv[k] += sLc->pinv[k];


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////










  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  cs_spfree(Bpinv);
  cs_spfree(W); cs_spfree(tW);
  cs_spfree(R); cs_spfree(Rinv); 
  if(lambda[0] == 0){
    cs_spfree(KRinv); cs_spfree(tWKRinv);
  }
  cs_spfree(tWKRinvW);
  // ttWKRinvW freed just after it is made so do not free here
  //free C just after making it/sLc - have room for a C in cs_remlcs_spfree(Ctmp);  
  cs_spfree(RHS); cs_spfree(tmpBLUXs); cs_spfree(BLUXs);

  cs_sfree(sLc);
  cs_nfree(Lc);

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
  delete [] newnu;
  delete [] dnu;
  delete [] trace;
  delete [] tugug;
  delete [] w;
  delete [] wchBd;
  delete [] con;
  delete [] cc;
  delete [] rfxlvls;


if(v[0] > 3){
  took = simple_toc(t);
  Rprintf("\n    %6.4f sec.: cpp post-REML freeing-up\n", took);
}

}
}

