#define _GREMLIN_H
#include "cs.h"
/* #include "R.h" included by cs.h */ 
#include "Rmath.h" 

#ifdef __cplusplus
extern "C" {
#endif




/* Expectation Maximization Algorithm:
     replaces elements in theta. Returns 1=success else 0 */
csi cs_em(const cs *BLUXs, const double *theta,
	csi nG, csi *rfxlvls, csi nb, csi *ndgeninv,
	cs **geninv, cs *Cinv);


cs *cs_chol2inv(const cs *L);

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
double cs_invR(const cs *C, const cs *A);
/* overwrites A with the inverse of the dense matrix C*/
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
/* FIXME not used cs *cs_initialize(double *x, int *p, int *i, int n, int m, int nzmax); */
/* allocate and fill a cs sparse matrix */
cs *cs_schur(const cs *A,  int split, const cs *beta);
/* forms the Schur complement for the dense mxm matrix: A_22-A_21%*%solve(A_11)%*%A_12 where submatrices are defined by split. Also overwrites beta_rr with A_21%*%solve(A_11) */

#ifdef __cplusplus
}
#endif

