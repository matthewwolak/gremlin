#include "gremlin.h"



/*
 note the trace of a product is equal to the sum of the element-by-element product
TODO FIXME check what to do if no ginverse associated with a parameter?!
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
XXX eqn. 2.44 is the score/gradient! for a varcomp

 current indexing faster than creating logical index matrix to multiply through Cinv for the purposes of subsetting/indexing
*/

/* theta overwritten with solution */
csi cs_em(const cs *BLUXs, double *theta, csi nG, csi *rfxlvls, csi nb){

  cs      *tmp_sln;
  double  *Bx;
  csi     g, k, si, qi, ei;
  if (!CS_CSC (BLUXs) || !nb || !theta) return (0);    // check arguments
  Bx = BLUXs->x;
 
  si = nb;
  for(g = 0; g < nG; g++){
    qi = rfxlvls[g];
    ei = si - 1 + qi;
//    tmp_sln = spalloc(qi, 1, qi, true, false);
//    for(k = 0; k < qi; k++){
//      tmp_sln
  }

  return(1);
// use  cs_idone() to finish and return?
}

