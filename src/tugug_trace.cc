#include "gremlin.h"

//#include "gremlincc.h"  // lists `gremlin.h` within it
// if NOT using clock functions simple_tic etc. then just include `gremlin.h`




// `tugugFun` function first, followed by `traceFun` function




/* return 1 if successful else returns 0
       tugug overwritten with output
       					 
   `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
*/
csi tugugFun(double *tugug, double *w, csi nG, csi *rfxlvls, csi *con,
        csi nb, csi *ndgeninv, cs **geninv, const cs *BLUXs
){

  double  *Bx;
  csi     g, i, k, si, qi, ei;

// for printing timings within
//double	t[2], took;
//simple_tic(t);

  if(!CS_CSC (BLUXs) || !nb) return (0);    // check arguments

  Bx = BLUXs->x;
  si = nb;

  //// `g` is the gth component of the G-structure to model
  //// `geninv` is the generalized inverse
  ////// (not the inverse of the G-matrix/(co)variance components)
  for(g = 0; g < nG; g++){
    tugug[g] = 0.0;
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    // skip parameter if it is fixed in the model (must do after calculate ei above)
    if(con[g] == 0){
      si = ei + 1;  // need to advance si first
      continue;  
    }

//TODO XXX for using covariance matrices, see Johnson & Thompson 1995 eqn 11a
    // Johnson & Thompson 1995 equations 9 & 10
    //// BUT t(BLUXs) %*% geninv == transpose of vector of geninv %*% BLUXs
    // `tugug` is the numerator of the 3rd term in the equation
    //// Above post-multiplied by BLUXs (... %*% BLUXs)
    if(ndgeninv[g] > 0){
      // Non-diagonal generalized inverses
      // crossproduct of BLUXs[si:ei] with geninv (t(BLUXs) %*% geninv %*% BLUXs)
      if(!CS_CSC(geninv[g])) error("geninv[%i] not CSC matrix\n", g);
      for(i = 0; i < geninv[g]->n; i++) w[i] = 0.0;  // clear w
      // Johnson & Thompson 1995 eqn 9a
      for(k = 0; k < geninv[g]->n; k++){
        for(i = geninv[g]->p[k]; i < geninv[g]->p[k+1]; i++){
          w[k] += geninv[g]->x[i] * Bx[geninv[g]->i[i]+si];
        }  // end for i
        tugug[g] += w[k] * Bx[k+si];
      }  // end for k
    } else{
        for(k = si; k <= ei; k++) tugug[g] += Bx[k] * Bx[k];  
      }  // end if/else geninv is non-diagonal
    si = ei + 1;

// print timing
/*
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: gradFun calculate tugug for %i of %i term", took, g, nG-1);
simple_tic(t);
*/

  }  // end for g in nG

  return(1);
}






////////////////////////////////////////////////////////////////////////////////




/*
 note the trace of a product is equal to the sum of the element-by-element product
  don't have to make `Cinv` for diagonal generalized inverses, just Cinv diagonals
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
*/

/* return 1 if successful else returns 0
       trace overwritten with output
*/
csi traceFun(double *trace, double *w,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv, cs **geninv,
	csi nsln, const cs *Lc, const csi *Pinv
){

  double  *Lx;
  csi     g, i, j, k, si, qi, ei, *Lp, *Li, minginviPerm;

// for printing timings within
/*
double	t[2], took;
simple_tic(t);
*/

  if(!nb) return (0);    // check arguments

  Lp = Lc->p; Li = Lc->i; Lx = Lc->x;
  si = nb;

  // `trace` = trace(Cinv_gg %*% geninv_gg])
  //// `g` is the gth component of the G-structure to model
  //// `geninv` is the generalized inverse
  ////// (not the inverse of the G-matrix/(co)variance components)
  for(g = 0; g < nG; g++){
    trace[g] = 0.0;
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    // DO NOT skip parameter if it is fixed in the model
    //// trace of all components is necessary to calculate gradient of residual

    // trace(geninv[g] %*% Cinv[si:ei, si:ei]) ...
    //// geninv[g] %*% Cinv[si:ei, si:ei]
    ////// each k column of Lc, create Cinv[si:ei, si:ei] elements
    //////// forward solve with Identity[, k] then backsolve
    for(k = si; k <= ei; k++){

// timing
/*
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  // just setting start time, no need to print here
  simple_tic(t);
}
*/

      // use w as working vector: first create Identity[, k]
      for(i = 0; i < nsln; i++) w[i] = 0.0;  // clear w

// print timing
/*
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: clear w at %i for %i of %i term", took,k,g,nG-1);
  simple_tic(t);
}
*/

      //// Lc is permuted, use t(P) %*% I[, k]
      w[Pinv[k]] += 1.0;              /* essentially `cs_ipvec` */
      // gr_cs_lltsolve(Lc, w, Pinv[k]); /* x = L\x then x = L'\x  */
      // forward solve (e.g., cs_lsolve)
      for(j = Pinv[k]; j < nsln; j++){
        if(w[j] != 0.0){
          w[j] /= Lx[Lp[j]];  // set diagonal (1 / L[k,k])
          // for loop to determine off-diagonal contributions
          for(i = Lp[j] + 1; i < Lp[j + 1]; i++){
            w[Li[i]] -= Lx[i] * w[j];
          }
        }  // end if NOT 0
      }

// print timing AND CHECK
/*
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: lsolve step %i for %i of %i term",took,k,g,nG-1);
  simple_tic(t);
}
*/

       // if Diagonal generalized inverse associated
       if(ndgeninv[g] == 0){
         // only need to find diagonal element for Cinv_ii[k] to calculate trace
         for(j = nsln - 1; j >= Pinv[k]; j--){
           for(i = Lp[j] + 1; i < Lp[j + 1]; i++){
             w[j] -= Lx[i] * w[Li[i]];
           }
           w[j] /= Lx[Lp[j]];
         }  
       } else{

           // calculate minimum row of geninverse[g][, k] after permuting it
           minginviPerm = nsln;
           for(i = geninv[g]->p[k - si]; i < geninv[g]->p[k - si + 1]; i++){
             j = Pinv[geninv[g]->i[i] + si];
             if(j < minginviPerm) minginviPerm = j;
           }

           // if non-diagonal generalized inverse
           //// need entire "column" of Cinv
           ////// until worked up to first row of non-zero in geninv[g][, k]
           for(j = nsln - 1; j >= minginviPerm; j--){
             for(i = Lp[j] + 1; i < Lp[j + 1]; i++){
               w[j] -= Lx[i] * w[Li[i]];
             }
             w[j] /= Lx[Lp[j]];
           }  
         }  // end if/else non-diag geninv[g]


// print timing AND CHECK
/*
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: ltsolve step %i for %i of %i term",took,k,g,nG-1);
  simple_tic(t);
}
*/

      // Now w holds a permuted ordering of Cinv[, k]
      if(ndgeninv[g] == 0){
        // if a diagonal geninv, trace=sum(diag(Cinv[si:ei, si:ei]))
        //// trace (numerator of 2nd) term in the equation
        trace[g] += w[Pinv[k]];  // essentially `cs_pvec` on w
      } else{
        // Contribution to the trace: calculate diagonal of matrix product
        //(geninv[g]%*%Cinv[si:ei, si:ei])[k,k] = geninv[g][k,]%*%Cinv[si:ei,k]
        //// trace (numerator of 2nd) term in the equation
        for(i = geninv[g]->p[k - si]; i < geninv[g]->p[k - si + 1]; i++){
          trace[g] += geninv[g]->x[i] * w[Pinv[geninv[g]->i[i] + si]];
        }  // end for i
      }  // end if/else non-diagonal geninv trace calculation
    }  // end for k
    si = ei + 1;

// print timing
/*
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: calculate trace for %i of %i term", took, g, nG-1);
simple_tic(t);
*/


  }  // end for g in nG


  return(1);
}



