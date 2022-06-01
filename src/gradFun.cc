#include "gremlin.h"


// `cs_gradFun` function first, followed by `cs_gradFun_fd` function



/*
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp
*/


/* return 1 if successful else returns 0
       dLdnu overwritten with output
       					 */
csi cs_gradFun(double *nu, double *dLdnu, 
        double *tugug, double *trace, csi *con,
	csi n, csi nG, csi *rfxlvls, csi nb,
	double sigma2e,    // 1.0 if lambda=FALSE
	csi thetaR, double *r      // 0 if lambda=TRUE
){

  int     lambda, nminfrfx;
  csi     g, k;

  if(!nb || !nu) return (0);    // check arguments
  if(thetaR != 0 && fabs(sigma2e - 1.00) < DBL_EPSILON) lambda = 0; else lambda = 1;

  nminfrfx = n - nb;
    if(lambda == 0) for(g = 0; g < nG; g++) nminfrfx -= rfxlvls[g];

  // First derivatives (gradient/score)
  // Johnson and Thompson 1995 don't use -0.5
  //// because likelihood is -2 log likelihood
  //// see `-2` on left-hand side of Johnson & Thompson eqn 3
  // Johnson and Thompson 1995: base to Appendix 2 eqn B3 and eqn 9a and 10a
  for(g = 0; g < nG; g++){
    dLdnu[g] = (con[g] == 0) ? 0.0 : (rfxlvls[g] / nu[g]);
  }  // end for g

  if(lambda == 1){
    for(g = 0; g < nG; g++){

      // Johnson and Thompson 1995 Appendix 2 eqn B3 and eqn 9a and 10a
      dLdnu[g] -= (1 / (nu[g] * nu[g])) * (trace[g] + tugug[g] / sigma2e);
      dLdnu[g] *= -0.5;
    }  // end for g

  } else{  // else when NOT lambda scale
    // Johnson and Thompson 1995 eqn 9b
    tugug[thetaR] = 0.0;
      for(k = 0; k < n; k++) tugug[thetaR] += r[k] * r[k];  // crossprod(residual)
//FIXME change `[thetaR]` below to be number of residual (co)variances
    dLdnu[thetaR] = (nminfrfx / nu[thetaR]) - (tugug[thetaR] / (nu[thetaR] * nu[thetaR]));
    for(g = 0; g < nG; g++){
      dLdnu[thetaR] += (1 / nu[thetaR]) * (trace[g] / nu[g]);
      // Johnson and Thompson 1995 eqn 9a and 10a
      dLdnu[g] -= (1 / (nu[g] * nu[g])) * (trace[g] + tugug[g]);
      dLdnu[g] *= -0.5;
    }  // end for g
    dLdnu[thetaR] *= -0.5;
  }  // end when NOT lambda scale

  // set any gradients for fixed parameters = 0.0
  for(g = 0; g < nG; g++) if(con[g] == 0) dLdnu[g] = 0.0;
//FIXME change `[thetaR]` below to be number of residual (co)variances
  if(con[thetaR] == 0) dLdnu[thetaR] = 0.0;

  return(1);
}







////////////////////////////////////////////////////////////////////////////////





/* return 1 if successful else returns 0
       dLdnu overwritten with output
fd = 0; backwards finite differences
fd = 1; central (both backwards and forward finite differences)
fd = 2; forward finite differences       
      					 */
csi cs_gradFun_fd(double *nu, csi fd, double h,
	double *dLdnu, double lL, csi *con,
	csi n, csi *dimZWG, csi nG, csi p, double *y,
	cs *Bpinv, cs *W, cs *tW, csi *rfxlvls, double rfxlL,
	csi *ndgeninv, cs **geninv, cs *KRinv,
	cs *Ctmp, cs *RHS, cs *tmpBLUXs, cs *BLUXs,
	css *sLc,
	double tyRinvy, // lambda=TRUE same every iteration else 0.0 when FALSE
	int nminffx, // lambda=TRUE same every iteration else 0
	int *nnzGRs, int *dimGRs, int *iGRs, csi lmbda
){

  int     g, k, si, dimM;
  double  tyPy, logDetC, sigma2e, loglik, denomSC;
          
  cs      *R, *Rinv, *tWKRinv, *tWKRinvW, *ttWKRinvW;

  cs*     *G = new cs*[nG];
  cs*     *Ginv = new cs*[nG];
  cs* 	  *KGinv = new cs*[nG];

  csn     *Lc;

  double *r = new double[dimZWG[2]];
  double *fxL = new double[p];
  double *fxU = new double[p];

  if(!lL || !nu) return (0);    // check arguments
  
  sigma2e = (lmbda == 1) ? 0.0 : 1.0;

  denomSC = (fd == 1) ? 2.0 : 1.0;

  // seed upper and lower log-likelihood vectors with current model log-likelihood
  //// so if using either backward or forward then will be subtracting from this
  for(g = 0; g < p; g++){
    fxL[g] = lL;
    fxU[g] = lL;
  }
  
  
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
  if(lmbda == 1){
    tWKRinvW = cs_multiply(tW, W);
  } else{
      // if lambda=FALSE    
      cs_kroneckerIupdate(Rinv, dimZWG[2], KRinv); 
      // Components of Meyer 1989 eqn 2
      tWKRinv = cs_multiply(tW, KRinv);
      tWKRinvW = cs_multiply(tWKRinv, W);
    }
    
  // Now take transpose of transpose to correctly order (don't ask why)
  ttWKRinvW = cs_transpose(tWKRinvW, true);
  cs_spfree(tWKRinvW);
  tWKRinvW = cs_transpose(ttWKRinvW, true);
  cs_spfree(ttWKRinvW);


  // Initial setup of G matrices
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
    // initial setup of KGinv inside gradFun
    if(ndgeninv[g] == 0){   // Diagonal matrix kronecker product: Ginv %x% I
      KGinv[g] = cs_kroneckerI(Ginv[g], dimZWG[4+g]);
    } else{  // generalized inverse kronecker product: Ginv %x% geninv
        KGinv[g] = cs_kroneckerA(Ginv[g], geninv[g]);
      }   
  }  // end for g





  //////////////////////////////////////////////////////////////////////////////
  // Skip most of below if there are no variance components other than residual
  if(nG > 0){
    // `g` is the gth component of the G-structure model
    for(g = 0; g < nG; g++){
      if(con[g] == 0){
        dLdnu[g] = 0.0;
      } else{
          // replace elements in G with finite difference pertubation(s)
          //// if either FORWARD or CENTRAL difference method
          if(fd > 0){
            cs_spfree(Ginv[g]);    
            G[g]->x[0] = nu[g] + h;   // TODO fix x[0] when G can be a matrix
            Ginv[g] = cs_inv(G[g]);
            Lc = cs_reml(n, dimZWG, nG, p, y,
              Bpinv, W, tW, rfxlvls, rfxlL,
              R, Rinv, G, Ginv, ndgeninv, geninv,
              KRinv, KGinv, tWKRinv, tWKRinvW, Ctmp,
              RHS, tmpBLUXs, BLUXs, r,
              sLc, 
              &tyPy, &logDetC, &sigma2e,
              tyRinvy,
              nminffx,
              &loglik,
              1, 0, 0, lmbda);      
if(loglik == 0.0){
  error("\nUnsuccessful REML calculation: finite difference gradient function component %i",
    g);
}
            cs_nfree(Lc);  
            fxL[g] = loglik;  // Either Forward or Central diff. lL
          }  // end if fd > 0  
      
          //// if either BACKWARD or CENTRAL difference method
          if(fd < 2){
            cs_spfree(Ginv[g]);    
            G[g]->x[0] = nu[g] - h;   // TODO fix x[0] when G can be a matrix
            Ginv[g] = cs_inv(G[g]);
            Lc = cs_reml(n, dimZWG, nG, p, y,
              Bpinv, W, tW, rfxlvls, rfxlL,
              R, Rinv, G, Ginv, ndgeninv, geninv,
              KRinv, KGinv, tWKRinv, tWKRinvW, Ctmp,
              RHS, tmpBLUXs, BLUXs, r,
              sLc, 
              &tyPy, &logDetC, &sigma2e,
              tyRinvy,
              nminffx,
              &loglik,
              1, 0, 0, lmbda);      
if(loglik == 0.0){
  error("\nUnsuccessful REML calculation: finite difference gradient function component %i",
    g);
}
            cs_nfree(Lc);      
            fxU[g] = loglik;  // Either Backward or Central diff. lL
          }  // end if fd < 2 
      
          // reset G
          G[g]->x[0] = nu[g];  // TODO fix x[0] when G can be a matrix

        }  // end when gth parameter NOT constrained

    }  // end for g

  }  // end if no G-structure variance components




  // Residual (co)variances when not on Lambda scale
  if(lmbda == 0){
//FIXME change `[p]` below to be number of residual (co)variances
    // `g` is the gth component of the R-structure model
    for(g = nG; g < p; g++){
      if(con[g] == 0){
        dLdnu[g] = 0.0;
      } else{
         // replace elements in R with finite difference pertubation(s)
          //// if either FORWARD or CENTRAL difference method
          if(fd > 0){
            cs_spfree(Rinv);
            cs_spfree(tWKRinv);
            cs_spfree(tWKRinvW);
                
            R->x[0] = nu[g] + h;   // TODO fix x[0] when R can be a matrix
            Rinv = cs_inv(R);
            cs_kroneckerIupdate(Rinv, dimZWG[2], KRinv); 
            // Components of Meyer 1989 eqn 2
            tWKRinv = cs_multiply(tW, KRinv);
            tWKRinvW = cs_multiply(tWKRinv, W);
            // Now take transpose of transpose to correctly order (don't ask why)
            ttWKRinvW = cs_transpose(tWKRinvW, true);
            cs_spfree(tWKRinvW);
            tWKRinvW = cs_transpose(ttWKRinvW, true);
            cs_spfree(ttWKRinvW);

            Lc = cs_reml(n, dimZWG, nG, p, y,
              Bpinv, W, tW, rfxlvls, rfxlL,
              R, Rinv, G, Ginv, ndgeninv, geninv,
              KRinv, KGinv, tWKRinv, tWKRinvW, Ctmp,
              RHS, tmpBLUXs, BLUXs, r,
              sLc, 
              &tyPy, &logDetC, &sigma2e,
              tyRinvy,
              nminffx,
              &loglik,
              1, 0, 0, lmbda);      
if(loglik == 0.0){
  error("\nUnsuccessful REML calculation: finite difference gradient function component %i",
    g);
}
            cs_nfree(Lc);      
            fxL[g] = loglik;  // Either Forward or Central diff. lL
          }  // end if fd > 0  
      
          //// if either BACKWARD or CENTRAL difference method
          if(fd < 2){
            cs_spfree(Rinv);
            cs_spfree(tWKRinv);
            cs_spfree(tWKRinvW);
                
            R->x[0] = nu[g] - h;   // TODO fix x[0] when R can be a matrix
            Rinv = cs_inv(R);
            cs_kroneckerIupdate(Rinv, dimZWG[2], KRinv); 
            // Components of Meyer 1989 eqn 2
            tWKRinv = cs_multiply(tW, KRinv);
            tWKRinvW = cs_multiply(tWKRinv, W);
            // Now take transpose of transpose to correctly order (don't ask why)
            ttWKRinvW = cs_transpose(tWKRinvW, true);
            cs_spfree(tWKRinvW);
            tWKRinvW = cs_transpose(ttWKRinvW, true);
            cs_spfree(ttWKRinvW);

            Lc = cs_reml(n, dimZWG, nG, p, y,
              Bpinv, W, tW, rfxlvls, rfxlL,
              R, Rinv, G, Ginv, ndgeninv, geninv,
              KRinv, KGinv, tWKRinv, tWKRinvW, Ctmp,
              RHS, tmpBLUXs, BLUXs, r,
              sLc, 
              &tyPy, &logDetC, &sigma2e,
              tyRinvy,
              nminffx,
              &loglik,
              1, 0, 0, lmbda);      
if(loglik == 0.0){
  error("\nUnsuccessful REML calculation: finite difference gradient function component %i",
    g);
}
            cs_nfree(Lc);      
            fxU[g] = loglik;  // Either Backward or Central diff. lL
          }  // end if fd < 2 
      
          // reset R
          R->x[0] = nu[g];  // TODO fix x[0] when R can be a matrix
      }  // end when parameter NOT constrained
      
    }  // end for g through R-structure components 

  }  // end when lambda FALSE


  // First derivatives (gradient/score)
  //// backward = [f(x) - f(x-h)] / h
  //// central = [f(x+h) - f(x-h)] / 2h
  //// forward = [f(x+h) - f(x)] / h
  // unpack/calculate derivativs of log-likelihood from differences 
  for(g = 0; g < p; g++){
    dLdnu[g] = (fxL[g] - fxU[g]) / (denomSC * h);  
  }
 
  //////////////////////////////////////////////////////////////////////////////
  // Cleanup:
  cs_spfree(R);
  cs_spfree(Rinv);
  if(lmbda == 0) cs_spfree(tWKRinv);
  cs_spfree(tWKRinvW);
  // Lc always "freed" just after making with cs_reml
  
  //
  for(g = 0; g < nG; g++){
    cs_spfree(G[g]);
    cs_spfree(Ginv[g]);
    cs_spfree(KGinv[g]);
  }
  //
  delete [] G;
  delete [] Ginv;
  delete [] KGinv;
  //  
  delete [] r;
  delete [] fxL;
  delete [] fxU;

 return(1);
}


