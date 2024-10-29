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
      if(!CS_CSC(geninv[g])) Rf_error("geninv[%i] not CSC matrix\n", g);
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
 NOTE the trace of a product is equal to the sum of the element-by-element
 product. Alternatively, think about matrix multiplication to only compute
 diagonals - this ends up being a dot product.
 
 Don't have to make `Cinv` for diagonal generalized inverses, just sum(Cinv)
 diagonals
XXX TODO see Knight 2008 thesis eqns 2.36 & 2.42 (and intermediates) for more generalized form of what is in Mrode (i.e., multivariate/covariance matrices instead of single varcomps)
*/

/* return 1 if successful else returns 0
       trace overwritten with output
*/

csi traceFun(double *trace,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv, cs **geninv,
	csi nsln, const cs *Cinv, const csi *Pinv, double *Cinv_ii
){

  csi     g, j, k, r, t, si, qi, ei, nz = 0, *Cp, *Ci, *P;
  cs      *Cinv_gg, *tCinv_gg;
  
// for printing timings within
/*
double	t[2], took;
simple_tic(t);
*/

  if(!nb) return (0);    // check arguments

  Cp = Cinv->p; Ci = Cinv->i;
  P = cs_pinv(Pinv, nsln);
  si = nb;

  // `trace` = trace(Cinv_gg %*% geninv_gg])
  //// `g` is the gth component of the G-structure to model
  //// `geninv` is the generalized inverse
  ////// (not the inverse of the G-matrix/(co)variance components)
  ////// trace (numerator of 2nd) term in the equation
  for(g = 0; g < nG; g++){
    trace[g] = 0.0;
    qi = rfxlvls[g];
    ei = si - 1 + qi;
    // DO NOT skip parameter if it is fixed in the model
    //// trace of all components is necessary to calculate gradient of residual


    if(ndgeninv[g] == 0){
      // if a diagonal geninv:
      //// trace=sum(diag(Cinv[si:ei, si:ei]))
      ////      =sum(Cinv_ii[si:ei])
      ////// Cinv_ii should be in original order (not order of Cholesky factor)
      for(k = si; k <= ei; k++) trace[g] += Cinv_ii[k];
       
      
    } else{
      // trace(geninv[g] %*% Cinv[si:ei, si:ei]) ...
      // Contribution to trace: sum diagonal of matrix (element-wise) product
      //(geninv[g]%*%Cinv[si:ei, si:ei])[k,k] = geninv[g][, k] * Cinv[si:ei, k]
      //// Use partial C-inverse and re-order so in same order as geninv
      ////// call this Cinv_gg     
      ////// "steal"/strip down much of cs_permute to subset and permute Cinv


      // calculate conservative number of non-zeroes in Cinv_gg
      //// actually gives number of non-zeroes for all rows: Cinv[0:nsln, j]
      for(k = si; k <=ei; k++){
        j = Pinv[k];
        nz += Cp[j+1] - Cp[j];
      }
      Cinv_gg = cs_spalloc(qi, qi, nz, true, false);

      
      // find kth column that should be in Cinv_gg, add it, and re-order rows
      nz = 0;
      for(k = 0; k < qi; k++){
        Cinv_gg->p[k] = nz;
        
        j = Pinv[ k + si ];  // in general column j is column Pinv[k] of Cinv
        for(t = Cp[j]; t < Cp[j+1]; t++){
          r = P[ Ci[t] ];    // in general row r of column j is row P[t] of Cinv
          // only add row r if it s between rows si and ei of re-ordered Cinv
          //// if r not between si and ei it is above/below sub-matrix Cinv_gg
          if((r >= si) && (r <= ei)){
            Cinv_gg->x[nz] = Cinv->x[t];
            Cinv_gg->i[nz++] = r - si;
          }  // end if
        }  // end for t
        
      }  // end for k
      Cinv_gg->p[qi] = nz;    // add last column count of Cinv_gg
      // order Cinv_gg (mostly make sure rows are sorted)
      tCinv_gg = cs_transpose(Cinv_gg, 1);
        cs_spfree(Cinv_gg);
      Cinv_gg = cs_transpose(tCinv_gg, 1);
        cs_spfree(tCinv_gg);



      // Element-wise/Hadamard/Dot product: accumulating into trace
      for(k = 0; k < qi; k++){
        j = geninv[g]->p[k];
        r = Cinv_gg->p[k];
        
        while((j < geninv[g]->p[k+1]) && (r < Cinv_gg->p[k+1])){
          t = geninv[g]->i[j];
          if(t == Cinv_gg->i[r]){
            trace[g] += geninv[g]->x[j++] * Cinv_gg->x[r++];
          } else (t > Cinv_gg->i[r]) ? r++ : j++ ;
        }  // end while
        
      }  // end for k
           
      cs_spfree(Cinv_gg);  // cleanup
      
    }  // end if/else non-diagonal ginverse
    si = ei + 1;

// print timing AND CHECK
/*
if((k == si) | (k ==ei-1)){
  took = simple_toc(t);
  Rprintf("\n\t\t\t    %6.6f sec.: ltsolve step %i for %i of %i term",took,k,g,nG-1);
  simple_tic(t);
}
*/

// print timing
/*
took = simple_toc(t);
Rprintf("\n\t\t    %6.6f sec.: calculate trace for %i of %i term", took, g, nG-1);
simple_tic(t);
*/


  }  // end for g in nG

  cs_free(P);
  
  return(1);
}


