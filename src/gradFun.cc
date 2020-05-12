#include "gremlincc.h"  // lists `gremlin.h` within it
// if NOT using clock functions simple_tic etc. then just include `gremlin.h`



/*
 note the trace of a product is equal to the sum of the element-by-element product
Consequently, don't have to make `Cinv`, just diagonals
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

*/

/* return 1 if successful else returns 0
       dLdnu overwritten with output
       Cinv_ii overwritten with diag(Cinv) */
csi cs_gradFun(double *nu, double *dLdnu, double *Cinv_ii,
	csi n, csi p, csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv,
	const cs *BLUXs, const cs *Lc, const csi *Pinv,
        double sigma2e,    // 1.0 if lambda=FALSE
	csi thetaR, double *r      // 0 if lambda=TRUE
){

  int     lambda, nminfrfx;
  double  *Bx;
  csi     g, i, j, k, nsln, si, qi, ei, qimax;

// for printing timings within gradFun
double	t[2], took;
simple_tic(t);


  if(!CS_CSC (BLUXs) || !nb || !nu) return (0);    // check arguments
  if(thetaR != 0 && fabs(sigma2e - 1.00) < DBL_EPSILON) lambda = 0; else lambda = 1;

  nsln = BLUXs->m;
  Bx = BLUXs->x;
  si = nb;
  nminfrfx = n - nb;
    if(lambda == 0) for(g = 0; g < nG; g++) nminfrfx -= rfxlvls[g];
  // `Cinv_ii` is the diagonals only of the C-inverse matrix
  // `trace` = trace(Cinv_gg %*% geninv_gg])
  // `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
  //// `g` is the gth component of the G-structure to model
  //// `geninv` is the generalized inverse
  ////// (not the inverse of the G-matrix/(co)variance components)
  double  *trace = new double[nG];
  if(lambda == 0) g = nG + 1; else g = nG;
  double  *tugug = new double[g];  // includes crossprod(residual) when !lambda
  double  *w = new double[nsln];

  // only need sln_g for non-diagonal generalized inverses, so make as large as that
  //// if no non-diagonal generalized inverse this makes sln_g size 1
  qimax = 1;  
  for(g = 0; g < nG; g++){
    if((ndgeninv[g] > 0) & (rfxlvls[g] > qimax)) qimax = rfxlvls[g];
  }
  double  *sln_g = new double[qimax];

  for(g = 0; g < nG; g++){
    trace[g] = 0.0;
    tugug[g] = 0.0;
    qi = rfxlvls[g];
    ei = si - 1 + qi;

//TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    // Johnson & Thompson 1995 equations 9 & 10
    //// BUT t(BLUXs) %*% geninv == transpose of vector of geninv %*% BLUXs
    for(i = 0; i < nsln; i++) w[i] = 0.0;  // clear w
    if(ndgeninv[g] > 0){
      // Non-diagonal generalized inverses
      // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv)
      if(!CS_CSC(geninv[g])) error("geninv[%i] not CSC matrix\n", g);
      for(i = 0; i < qi; i++) sln_g[i] = Bx[si + i];
      // Johnson & Thompson 1995 eqn 9a
      cs_gaxpy(geninv[g], sln_g, w);
      for(i = 0; i < qi; i++) sln_g[i] = 0.0;
    } else{
      for(i = 0; i < qi; i++) w[i] = Bx[si + i];  
    }  // end if/else geninv is non-diagonal

    // `tugug` is the numerator of the 3rd term in the equation
    //// Above (w) post-multiplied by BLUXs (... %*% BLUXs)
    for(k = 0; k < qi; k++) tugug[g] += w[k] * Bx[si + k];



// print timing
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: calculate tugug for %i of %i term", took, g, nG-1);
simple_tic(t);



    // ... + trace(geninv[g] %*% Cinv[si:ei, si:ei]) ...
    //// geninv[g] %*% Cinv[si:ei, si:ei]
    ////// each k column of Lc, create Cinv[si:ei, si:ei] elements
    //////// forward solve with Identity[, k] then backsolve
    for(k = si; k <= ei; k++){
// timing
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  // just setting start time, no need to print here
  simple_tic(t);
}
    // use w as working vector: first create Identity[, k]
      for(i = 0; i < nsln; i++) w[i] = 0.0;  // clear w
 // print timing
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: clear w at %i for %i of %i term", took,k,g,nG-1);
  simple_tic(t);
}
     //// Lc is permuted, use t(P) %*% I[, k]
      w[Pinv[k]] += 1.0;              /* essentially `cs_ipvec` */
      gr_cs_lltsolve(Lc, w, Pinv[k]); /* x = L\x then x = L'\x  */
// print timing AND CHECK
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: gr_cs_lltsolve %i for %i of %i term",took,k,g,nG-1);
  simple_tic(t);
}
      // Now w holds a permuted ordering of Cinv[, k]
      // write sampling variance for kth BLUP
      Cinv_ii[k] = w[Pinv[k]];	      /* essentially `cs_pvec` */





      if(ndgeninv[g] == 0){
        // if a diagonal geninv, trace=sum(diag(Cinv[si:ei, si:ei]))
        //// trace (numerator of 2nd) term in the equation
        trace[g] += Cinv_ii[k];
      } else{
        // Contribution to the trace: calculate diagonal of matrix product
        //(geninv[g]%*%Cinv[si:ei, si:ei])[k,k] = geninv[g][k,]%*%Cinv[si:ei,k]
        //// trace (numerator of 2nd) term in the equation
        for(i = geninv[g]->p[k-si]; i < geninv[g]->p[k-si+1]; i++){
          j = geninv[g]->i[i] + si;
          trace[g] += geninv[g]->x[i] * Cinv_ii[j];
        }  // end for i
      }  // end if/else non-diagonal geninv trace calculation
    }  // end for k
    si = ei + 1;


// print timing
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: calculate trace for %i of %i term", took, g, nG-1);
simple_tic(t);



  }  // end for g in nG






  // First derivatives (gradient/score)
  // Johnson and Thompson 1995 don't use -0.5
  //// because likelihood is -2 log likelihood
  //// see `-2` on left-hand side of Johnson & Thompson eqn 3
  // Johnson and Thompson 1995: base to Appendix 2 eqn B3 and eqn 9a and 10a
  for(g = 0; g < nG; g++){
    dLdnu[g] = (rfxlvls[g] / nu[g]);
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



  delete [] w;
  delete [] tugug;
  delete [] trace;
  delete [] sln_g;

  return(1);
}

















////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
 note the trace of a product is equal to the sum of the element-by-element product
Consequently, don't have to make `Cinv`, just diagonals
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

*/

/* return 1 if successful else returns 0
       dLdnu overwritten with output
       Cinv_ii overwritten with diag(Cinv) */
csi cs_gradFun2(double *nu, double *dLdnu, double *Cinv_ii,
	csi n, csi p, csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv,
	const cs *BLUXs, const cs *Lc, const csi *Pinv,
        double sigma2e,    // 1.0 if lambda=FALSE
	csi thetaR, double *r      // 0 if lambda=TRUE
){

  int     lambda, nminfrfx;
  double  *Bx;
  csi     g, i, j, k, nsln, si, qi, ei, qimax;

// for printing timings within gradFun
double	t[2], took;
simple_tic(t);


  if(!CS_CSC (BLUXs) || !nb || !nu) return (0);    // check arguments
  if(thetaR != 0 && fabs(sigma2e - 1.00) < DBL_EPSILON) lambda = 0; else lambda = 1;

  nsln = BLUXs->m;
  Bx = BLUXs->x;
  si = nb;
  nminfrfx = n - nb;
    if(lambda == 0) for(g = 0; g < nG; g++) nminfrfx -= rfxlvls[g];
  // `Cinv_ii` is the diagonals only of the C-inverse matrix
  // `trace` = trace(Cinv_gg %*% geninv_gg])
  // `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
  //// `g` is the gth component of the G-structure to model
  //// `geninv` is the generalized inverse
  ////// (not the inverse of the G-matrix/(co)variance components)
  double  *trace = new double[nG];
  if(lambda == 0) g = nG + 1; else g = nG;
  double  *tugug = new double[g];  // includes crossprod(residual) when !lambda
  double  *w = new double[nsln];

  // only need sln_g for non-diagonal generalized inverses, so make as large as that
  //// if no non-diagonal generalized inverse this makes sln_g size 1
  qimax = 1;  
  for(g = 0; g < nG; g++){
    if((ndgeninv[g] > 0) & (rfxlvls[g] > qimax)) qimax = rfxlvls[g];
  }
  double  *sln_g = new double[qimax];

  for(g = 0; g < nG; g++){
    trace[g] = 0.0;
    tugug[g] = 0.0;
    qi = rfxlvls[g];
    ei = si - 1 + qi;

//TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    // Johnson & Thompson 1995 equations 9 & 10
    //// BUT t(BLUXs) %*% geninv == transpose of vector of geninv %*% BLUXs
    for(i = 0; i < nsln; i++) w[i] = 0.0;  // clear w
    if(ndgeninv[g] > 0){
      // Non-diagonal generalized inverses
      // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv)
      if(!CS_CSC(geninv[g])) error("geninv[%i] not CSC matrix\n", g);
      for(i = 0; i < qi; i++) sln_g[i] = Bx[si + i];
      // Johnson & Thompson 1995 eqn 9a
      cs_gaxpy(geninv[g], sln_g, w);
      for(i = 0; i < qi; i++) sln_g[i] = 0.0;
    } else{
      for(i = 0; i < qi; i++) w[i] = Bx[si + i];  
    }  // end if/else geninv is non-diagonal

    // `tugug` is the numerator of the 3rd term in the equation
    //// Above (w) post-multiplied by BLUXs (... %*% BLUXs)
    for(k = 0; k < qi; k++) tugug[g] += w[k] * Bx[si + k];



// print timing
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: calculate tugug for %i of %i term", took, g, nG-1);
simple_tic(t);



    // ... + trace(geninv[g] %*% Cinv[si:ei, si:ei]) ...
    //// geninv[g] %*% Cinv[si:ei, si:ei]
    ////// each k column of Lc, create Cinv[si:ei, si:ei] elements
    //////// forward solve with Identity[, k] then backsolve
    for(k = si; k <= ei; k++){
// timing
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  // just setting start time, no need to print here
  simple_tic(t);
}
    // use w as working vector: first create Identity[, k]
      for(i = 0; i < nsln; i++) w[i] = 0.0;  // clear w
 // print timing
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: clear w at %i for %i of %i term", took,k,g,nG-1);
  simple_tic(t);
}
     //// Lc is permuted, use t(P) %*% I[, k]
      w[Pinv[k]] += 1.0;              /* essentially `cs_ipvec` */
      gr_cs_lltsolve(Lc, w, Pinv[k]); /* x = L\x then x = L'\x  */
// print timing AND CHECK
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: gr_cs_lltsolve %i for %i of %i term",took,k,g,nG-1);
  simple_tic(t);
}
      // Now w holds a permuted ordering of Cinv[, k]
      // write sampling variance for kth BLUP
      Cinv_ii[k] = w[Pinv[k]];	      /* essentially `cs_pvec` */





      if(ndgeninv[g] == 0){
        // if a diagonal geninv, trace=sum(diag(Cinv[si:ei, si:ei]))
        //// trace (numerator of 2nd) term in the equation
        trace[g] += Cinv_ii[k];
      } else{
        // Contribution to the trace: calculate diagonal of matrix product
        //(geninv[g]%*%Cinv[si:ei, si:ei])[k,k] = geninv[g][k,]%*%Cinv[si:ei,k]
        //// trace (numerator of 2nd) term in the equation
        for(i = geninv[g]->p[k-si]; i < geninv[g]->p[k-si+1]; i++){
          j = geninv[g]->i[i] + si;
          trace[g] += geninv[g]->x[i] * Cinv_ii[j];
        }  // end for i
      }  // end if/else non-diagonal geninv trace calculation
    }  // end for k
    si = ei + 1;


// print timing
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: calculate trace for %i of %i term", took, g, nG-1);
simple_tic(t);



  }  // end for g in nG






  // First derivatives (gradient/score)
  // Johnson and Thompson 1995 don't use -0.5
  //// because likelihood is -2 log likelihood
  //// see `-2` on left-hand side of Johnson & Thompson eqn 3
  // Johnson and Thompson 1995: base to Appendix 2 eqn B3 and eqn 9a and 10a
  for(g = 0; g < nG; g++){
    dLdnu[g] = (rfxlvls[g] / nu[g]);
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



  delete [] w;
  delete [] tugug;
  delete [] trace;
  delete [] sln_g;

  return(1);
}

