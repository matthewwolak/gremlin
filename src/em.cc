#include "gremlin.h"


// EM refs: Hofer 1998 eqn 10-12
// note Hofer eqn 12 has no sigma2e in last term of non-residual formula
//// ?mistake? in Mrode 2005 (p. 241-245), which does include sigma2e

/*
 note the trace of a product is equal to the sum of the element-by-element product
Consequently, don't have to make `Cinv`, just diagonals
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

*/

/* nu overwritten with solution */
/* Cinv_ii overwritten with diag(Cinv) */
csi cs_em(const cs *BLUXs, double *nu, double *Cinv_ii,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv, const cs *Lc, const csi *Pinv
){

  double  *Bx, r, o, tr;
  csi     g, i, j, k, n, si, qi, ei;
  if(!CS_CSC (BLUXs) || !nb || !nu) return (0);    // check arguments
  n = BLUXs->m;
  double  *tmp_sln = new double[n];
  double  *w = new double[n];
    for(k = 0; k < n; k++){
      tmp_sln[k] = 0.0;
      w[k] = 0.0;
    }
  Bx = BLUXs->x;
  si = nb;
  r = nu[nG];                       //   residual var

  for(g = 0; g < nG; g++){
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    o = 0.0; tr = 0.0;

    // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv)
    //// Diagonal geninv first
    if(ndgeninv[g] > 0){
      if(!CS_CSC(geninv[g])) error("geninv[%i] not CSC matrix\n", g);
      //// Non-diagonal generalized inverses
      ////// Column-by-column of geninv
      for(k = 0; k < qi; k++){
        j = k + si;                    // index in solution vector
        ////// Each row (i) in col k of geninv
        for(i = geninv[g]->p[k]; i < geninv[g]->p[k+1]; i++){
          tmp_sln[j] += Bx[geninv[g]->i[i] + si] * geninv[g]->x[i];
        }  // end for i (rows of geninv column k)
      }  // end for k (columns of geninv)
    } else{
        //// Diagonal 
        ////// fill in tmp_sln with BLUXs
        for(k = si; k <= ei; k++) tmp_sln[k] = BLUXs->x[k];
    }  // end if/else geninv is non-diagonal


    // Above post-multiplied by BLUXs (... %*% BLUXs)
    for(k = si; k <= ei; k++){
        o += tmp_sln[k] * Bx[k];
    }

    // ... + trace(geninv[g] %*% Cinv[si:ei, si:ei]) ...
    //// geninv[g] %*% Cinv[si:ei, si:ei]
      ////// each k column of Lc, create Cinv[si:ei, si:ei] elements
      //////// forward solve with Identity[, k] then backsolve
      for(k = si; k <= ei; k++){
        // use tmp_sln as working vector: first zero out and create Identity[, k]
        //// however, Lc is permuted, so need to use t(P) %*% I[, k]
        for(i = 0; i < n; i++) tmp_sln[i] = 0.0;
        tmp_sln[Pinv[k]] += 1.0;              /* essentially `cs_ipvec` */
        cs_lsolve (Lc, tmp_sln);                /* x = L\x */      
        cs_ltsolve(Lc, tmp_sln);                /* x = L'\x */
        // Now tmp_sln holds a permuted ordering of Cinv[, k]
//TODO combine pvec (take code out of function) and multiplication in next step
        cs_pvec(Pinv, tmp_sln, w, n);         /* b = P'*x */
        // w holds Cinv[, k]

        // write sampling variance for kth BLUP
        Cinv_ii[k] = w[k];

        if(ndgeninv[g] == 0){
          // if a diagonal geninv, trace=sum(diag(Cinv[si:ei, si:ei]))
          tr += w[k];
        } else{
          // Contribution to the trace: calculate diagonal of matrix product
          //(geninv[g]%*%Cinv[si:ei, si:ei])[k,k] = geninv[g][k,]%*%Cinv[si:ei,k]
          for(i = geninv[g]->p[k-si]; i < geninv[g]->p[k-si+1]; i++){
            j = geninv[g]->i[i] + si;
            tr += geninv[g]->x[i] * w[j];
          }  // end for i
        }  // end if/else non-diagonal geninv trace calculation
      }  // end for k


    // (first term + trace ) / qi
    nu[g] = (o + tr ) / qi;
    si = ei + 1;
    // reset tmp_sln to zero for next g
    if(g < (nG-1)) for(k = 0; k < n; k++) tmp_sln[k] = 0.0;
  }  // end for g in nG



  //  XXX	XXX  Residual Done outside of function!!  XXX	XXX



   delete [] tmp_sln;
   delete [] w;

  return(1);
}

