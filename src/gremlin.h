#define _GREMLIN_H

#include "cs.h"
/* #include "R.h" included by cs.h */ 
#include "Rmath.h" 

#ifdef __cplusplus
extern "C" {
#endif




/* returns log-likelihood if successful, 0 if not
   Replaces all objects that are updated based on changed `nu` values */
csn *cs_reml(csi n, csi *dimZWG, csi nG, csi p, double *y,
	cs *Bpinv, cs *W, cs *tW, csi *rfxlvls, double rfxlL,
	cs *R, cs *Rinv, cs **G, cs  **Ginv, csi *ndgeninv, cs **geninv,
	cs *KRinv, cs **KGinv, cs *tWKRinv, cs *tWKRinvW, cs *Ctmp,
	cs *RHS, cs *tmpBLUXs, cs *BLUXs, double *res,
	css *sLc,
	double *tyPy, double *logDetC, double *sigma2e, double *loglik,
	csi i, csi v, csi vitout);


/* Average Information Algorithm:
    returns replaces AI */
cs *cs_ai(cs *BLUXs, cs **Ginv,
        cs *R, cs *KRinv, cs *tWKRinv,
        double *rory,  // residuals if lambda=FALSE else y if lambda=TRUE
        cs *W, cs *tW, csi n, csi p, csi nG, csi *rfxlvls, csi nb,
	cs *Lc, csi *Pinv,
	csi thetaR,       // 0 if lambda=TRUE
        double sigma2e);    // 1.0 if lambda=FALSE


 
/* Analytical Gradient/Score (first derivative) function
     return 1 if successful else returns 0
     dLdnu overwritten with output */
csi cs_gradFun(double *nu, double *dLdnu, 
        double *tugug, double *trace, csi *con,
	csi n, csi nG, csi *rfxlvls, csi nb,
	double sigma2e,    // 1.0 if lambda=FALSE
	csi thetaR, double *r);      // 0 if lambda=TRUE


/* Finite difference Gradient/Score (first derivative) function
     return 1 if successful else returns 0
     dLdnu overwritten with output */
csi cs_gradFun_fd(double *nu, csi fd, double *dLdnu, double lL, csi *con,
	csi n, csi *dimZWG, csi nG, csi p, double *y,
	cs *Bpinv, cs *W, cs *tW, csi *rfxlvls, double rfxlL,
	csi *ndgeninv, cs **geninv, cs *KRinv,
	cs *Ctmp, cs *RHS, cs *tmpBLUXs, cs *BLUXs,
	css *sLc, 
	int *nnzGRs, int *dimGRs, int *iGRs);

/* solve Ax=k where Lx=b, L'b=k, and x, b, and k are dense.
   x=b on input, solution on output. */
// combines cs_lsolve followed by cs_ltsolve
//// start at element k of x
csi gr_cs_lltsolve (const cs *L, double *x, csi k);


// replaces elements in Cinv_ii inverse diagonals as double
//// Returns 1=success else 0
csi cs_chol2inv_ii(const cs *L, const csi *Pinv, double *Cinv_ii, int r);  


/* return 1 if successful else returns 0
       tugug overwritten with output
       					 
   `tugug` = t(u_gg) %*% geninv_gg %*% u_gg
*/
csi tugugFun(double *tugug, double *w, csi nG, csi *rfxlvls, csi *con,
        csi nb, csi *ndgeninv, cs **geninv, const cs *BLUXs);

/* return 1 if successful else returns 0
       trace overwritten with output
*/
csi traceFun(double *trace, double *w,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv, cs **geninv,
	csi nsln, const cs *Lc, const csi *Pinv);

/* B  returned = the reduced matrix of A according to drop */
cs *cs_droprowcol(const cs *A, csi *drop);


/* 1 returned if successful, else NULL
   `newnu` replaced with next set of nu parameters   */
csi qNewtRhap(double *nu, double *newnu, double *dLdnu, const cs *A,
	      csi p, csi *con, csi *wchBd, double f, double *ezero, csi v);


/*******************************************************************/
/* Below are functions from MCMCglmm-2.25 by Jarrod Hadfield       */
/*******************************************************************/
cs *cs_cbind(const cs *A, const cs *B);
/* Returns the two matrices A and B column bound*/
void cs_cov2cor(const cs *A);
/* transforms a dense covariance matrix A into a correlation matrix  */
cs *cs_directsum(cs **KGinv, int nG, int nGR);
/* returns the direct sum of KGinv[1] KGinv[2] ... KGinv[nG]*/
void cs_directsumupdate(cs **KGinv, int nG, int nGR, const cs *C);
/* overwrites C with the direct sum of KGinv[1] KGinv[2] ... KGinv[nG]*/
cs *cs_inv(const cs *C);
/* returns the inverse of the dense matrix C*/
cs *cs_kroneckerA(const cs *G, const cs *A);
/* forms the kronecker product of G and A*/
void cs_kroneckerAupdate(const cs *G, const cs *A, const cs *C);
/* overwrites C with the kronecker product of G and A*/
cs *cs_kroneckerI(const cs *A, int nI);
/* forms the kronecker product of the dense matrix A and an identity matrix with dimension nI*/
void cs_kroneckerIupdate(const cs *A, int nI, const cs*C);
/* overwrites C with the kronecker product of the dense matrix A and an identity matrix with dimension nI*/
cs *cs_omega(cs **KGinv, int nG, cs *pvB);
/* returns the direct sum of pvB KGinv[1] KGinv[2] ... KGinv[nG] */
void cs_omegaupdate(cs **KGinv, int nG, cs *pvB, const cs *C);
/* overwrites C with the direct sum of pvB KGinv[1] KGinv[2] ... KGinv[nG] */
cs *cs_schur(const cs *A,  int split, const cs *beta);
/* forms the Schur complement for the dense mxm matrix: A_22-A_21%*%solve(A_11)%*%A_12 where submatrices are defined by split. Also overwrites beta_rr with A_21%*%solve(A_11) */

#ifdef __cplusplus
}
#endif

