#include "gremlin.h"



/*
 note the trace of a product is equal to the sum of the element-by-element product
TODO FIXME check what to do if no ginverse associated with a parameter?!
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

 current indexing faster than creating logical index matrix to multiply through Cinv for the purposes of subsetting/indexing
*/

/* theta overwritten with solution */
csi cs_em(const cs *BLUXs, const double *theta,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv, cs *Cinv
){

  double  *Bx, *tmp_sln, roqi, o;
  cs      *preTR, *preTRprod;
  csi     g, i, j, k, m, si, qi, ei;
  if(!CS_CSC (BLUXs) || !nb || !theta) return (0);    // check arguments
  Bx = tmp_sln = BLUXs->x;

  // Find location of rfxlvls with maximum number of levels
//FIXME: DELETE DELETE if don't end up using!!!
  m = 0;
  if(nG > 0){
    if(nG > 1){
      for(g = 1; g < nG; g++) m = (rfxlvls[g] < rfxlvls[m]) ? m : g;
    }
  }
//XXX  DELETE??   XXX

  si = nb; o = 0.0;
  for(g = 0; g < nG; g++){
    qi = rfxlvls[g];
    roqi = theta[nG];                                 //   residual var / qi
    ei = si - 1 + qi;


    // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv)
    //// Diagonal geninv first
    if(ndgeninv[g] == 0){
      for(k = si; k <= ei; k++){
        tmp_sln[k] *= qi;
      }
    } else{
      //// Non-diagonal generalized inverses
      ////// Zero gth range
      for(k = si; k <= ei; k++) tmp_sln[k] = 0.0;
      ////// Column-by-column of geninv
      for(k = 0; k < qi; k++){
        ////// Each row (i) in col k of geninv
        for(i = geninv[g]->p[k]; i < geninv[g]->p[k+1]; i++){
          ////// Each element in Bx
          for(j = si; j <= ei; j++){
            tmp_sln[j] += Bx[j] * geninv[g]->x[i];
          }  // end for j (elements in BLUXs)
        }  // end for i (rows of geninv column k)
      }  // end for k (columns of geninv)
    }  // end if/else whether geninv is diagonal

    // Above post-multiplied by BLUXs (... %*% BLUXs)
    for(k = si; k <= ei; k++){
      for(i = si; i <= ei; i++){
        o += tmp_sln[k] * Bx[i];
      }  // end for i (rows of BLUXs)
    }  // end for k (cols of tmp_sln)


    // ... + tr(geninv[g] %*% Cinv[si:ei, si:ei]) * residual variance / qi
    preTR = cs_spalloc(qi, qi, Cinv->p[ei+1] - Cinv->p[si], true, false);
      i = 0;
      for(k = Cinv->p[si]; k < Cinv->p[ei+1]; k++){
        preTR->i[i] = Cinv->i[k];
        preTR->x[i] = Cinv->x[k];
        i++;
      }
      for(k = 0; k <= qi; k++){
        preTR->p[k] = Cinv->p[si+k];
      }
    preTRprod = cs_multiply(geninv[g], preTR);// * roqi

    si = ei + 1;
  }  // end for g in nG

cs_spfree(preTR); cs_spfree(preTRprod);

  return(1);
// use  cs_idone() to finish and return?
}

