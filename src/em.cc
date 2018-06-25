#include "gremlin.h"



/*
 note the trace of a product is equal to the sum of the element-by-element product
TODO FIXME check what to do if no ginverse associated with a parameter?!
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

*/

/* theta overwritten with solution */
csi cs_em(const cs *BLUXs, double *theta,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv, cs *Cinv
){

  double  *Bx, r, o, tr;
  cs      *preTR, *preTRprod;
  csi     g, i, j, k, m, si, qi, ei;
  if(!CS_CSC (BLUXs) || !nb || !theta) return (0);    // check arguments
  double  *tmp_sln = new double[BLUXs->m];
    for(k = 0; k < BLUXs->m; k++) tmp_sln[k] = 0.0;
  Bx = BLUXs->x;
  si = nb;
  r = theta[nG];                       //   residual var

  for(g = 0; g < nG; g++){
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    o = 0.0; tr = 0.0;

    // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv)
    //// Diagonal geninv first
    if(ndgeninv[g] > 0){
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


    // ... + trace(geninv[g] %*% Cinv[si:ei, si:ei]) * residual variance / qi
    preTR = cs_spalloc(qi, qi, qi*qi, true, false);
      i = 0;
      preTR->p[0] = 0;
      for(k = si; k <= ei; k++){
        for(j = Cinv->p[k]; j < Cinv->p[k+1]; j++){
          if( (Cinv->i[j] >= si) && (Cinv->i[j] <= ei) ){
            preTR->i[i] = Cinv->i[j] - si;
            preTR->x[i] = Cinv->x[j];
            i++;
          }  // end if
        }  // end for j
        preTR->p[k-si+1] = i;
      }  // end for k
      preTR->nzmax = i;

      if(ndgeninv[g] == 0){
        preTRprod = cs_spalloc(preTR->m, preTR->n, preTR->nzmax, true, false);
          preTRprod->i = preTR->i; preTRprod->p = preTR->p;
          preTRprod->x = preTR->x;
      } else{
        if(!CS_CSC(geninv[g])) error("geninv[%i] not CSC matrix\n", g);
        preTRprod = cs_multiply(geninv[g], preTR);
      }
//      cs_spfree(preTR);
//      preTR = cs_transpose(preTRprod, 1);
//      cs_spfree(preTRprod);
//      preTRprod = cs_transpose(preTR, 1);

      // now tr(preTRprod)
      for(k = 0; k < qi; k++){
        for(j = preTRprod->p[k]; j < preTRprod->p[k+1]; j++){
          if(preTRprod->i[j] == k){
            tr += preTRprod->x[j];
            break;
          }
        }  // end for j
      }  // end for k
//      cs_spfree(preTRprod);
//      cs_spfree(preTR);

    // (first term + trace * r ) / qi
    theta[g] = (o + tr * r ) / qi;
    si = ei + 1;
  }  // end for g in nG



  //  XXX	XXX  Residual Done outside of function!!  XXX	XXX



  delete [] tmp_sln;
  return(1);
// use  cs_idone() to finish and return?
}

