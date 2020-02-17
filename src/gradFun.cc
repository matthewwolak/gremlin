#include "gremlin.h"


/*
 note the trace of a product is equal to the sum of the element-by-element product
Consequently, don't have to make `Cinv`, just diagonals
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

*/

/* return 1 if successful else returns 0
/* dLdnu overwritten with output */
/* Cinv_ii overwritten with diag(Cinv) */
csi cs_gradFun(double *nu, double *dLdnu, double *Cinv_ii,
	csi n, csi p, csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv,
	const cs *BLUXs, const cs *Lc, const csi *Pinv,
        double sigma2e,    // 1.0 if lambda=FALSE
	csi thetaR, double *r,      // 0 if lambda=TRUE
	double ezero
){

  int     lambda, nminfrfx;
  double  *Bx;
  csi     g, i, j, k, nsln, si, qi, ei;

  if(!CS_CSC (BLUXs) || !nb || !nu) return (0);    // check arguments
  if(thetaR != 0 && fabs(sigma2e - 1.00) < ezero) lambda = 0; else lambda = 1;

  nsln = BLUXs->m;
  Bx = BLUXs->x;
  si = nb;
  nminfrfx = n - nb;
    if(lambda == 0) for(g = 0; g < nG; g++) nminfrfx -= rfxlvls[g];
  // `Cinv_ii` is the diagonals only of the C-inverse matrix
  // `trace` = trace[Cinv_gg %*% geninv_gg]
  // `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
  //// `g` is the gth component of the G-structure to model
  //// `geninv` is the generalized inverse
  ////// (not the inverse of the G-matrix/(co)variance components)
  double  *trace = new double[nG];
  double  *tugug = new double[nG];  // includes crossprod(residual) when !lambda
  double  *tmp_sln = new double[nsln];
  double  *w = new double[nsln];
    for(k = 0; k < n; k++){
      tmp_sln[k] = 0.0;
      w[k] = 0.0;
    }

  for(g = 0; g < nG; g++){
    trace[g] = 0.0;
    tugug[g] = 0.0;
    qi = rfxlvls[g];
    ei = si - 1 + qi;

//TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    // Johnson & Thompson 1995 equations 9 & 10
    if(ndgeninv[g] > 0){
      // Non-diagonal generalized inverses
      // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv)
      if(!CS_CSC(geninv[g])) error("geninv[%i] not CSC matrix\n", g);
      // Johnson & Thompson 1995 eqn 9a
      ////// Column-by-column of geninv
      for(k = 0; k < qi; k++){
        j = k + si;                    // index in solution vector
        ////// Each row (i) in col k of geninv
        for(i = geninv[g]->p[k]; i < geninv[g]->p[k+1]; i++){
          tmp_sln[j] += Bx[geninv[g]->i[i] + si] * geninv[g]->x[i];
        }  // end for i (rows of geninv column k)
      }  // end for k (columns of geninv)

    } else{
        // Diagonal 
        //// fill in tmp_sln with BLUXs
        for(k = si; k <= ei; k++) tmp_sln[k] = Bx[k];
    }  // end if/else geninv is non-diagonal

    // `tugug` is the numerator of the 3rd term in the equation
    //// Above (tmp_sln) post-multiplied by BLUXs (... %*% BLUXs)
    for(k = si; k <= ei; k++){
        tugug[g] += tmp_sln[k] * Bx[k];
    }



    // ... + trace(geninv[g] %*% Cinv[si:ei, si:ei]) ...
    //// geninv[g] %*% Cinv[si:ei, si:ei]
    ////// each k column of Lc, create Cinv[si:ei, si:ei] elements
    //////// forward solve with Identity[, k] then backsolve
    for(k = si; k <= ei; k++){
    // use tmp_sln as working vector: first zero out and create Identity[, k]
    //// however, Lc is permuted, so need to use t(P) %*% I[, k]
      for(i = si; i < ei; i++) tmp_sln[i] = 0.0;  // clear tmp_sln
      tmp_sln[Pinv[k]] += 1.0;              /* essentially `cs_ipvec` */
      cs_lsolve (Lc, tmp_sln);              /* x = L\x */      
      cs_ltsolve(Lc, tmp_sln);              /* x = L'\x */
      // Now tmp_sln holds a permuted ordering of Cinv[, k]
//TODO combine pvec (take code out of function) and multiplication in next step
      cs_pvec(Pinv, tmp_sln, w, n);         /* b = P'*x */
      // w holds Cinv[, k]

      // write sampling variance for kth BLUP
      Cinv_ii[k] = w[k];

      if(ndgeninv[g] == 0){
        // if a diagonal geninv, trace=sum(diag(Cinv[si:ei, si:ei]))
        //// trace (numerator of 2nd) term in the equation
        trace[g] += w[k];
      } else{
        // Contribution to the trace: calculate diagonal of matrix product
        //(geninv[g]%*%Cinv[si:ei, si:ei])[k,k] = geninv[g][k,]%*%Cinv[si:ei,k]
        //// trace (numerator of 2nd) term in the equation
        for(i = geninv[g]->p[k-si]; i < geninv[g]->p[k-si+1]; i++){
          j = geninv[g]->i[i] + si;
          trace[g] += geninv[g]->x[i] * w[j];
        }  // end for i
      }  // end if/else non-diagonal geninv trace calculation
    }  // end for k

    si = ei + 1;
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
  delete [] tmp_sln;
  delete [] tugug;
  delete [] trace;

  return(1);
}

