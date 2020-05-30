#include "gremlin.h"

/* B  returned = the reduced matrix of A according to drop */
cs *cs_droprowcol(const cs *A, csi *drop)
{

  double  *Ax;
  cs      *B;
  csi     r, c, cnt, cd, rd, pr, *Ai, *Ap;

  if(!CS_CSC(A) || !drop) return(0);    // check arguments

  Ax = A->x; Ai = A->i; Ap = A->p;

  pr = A->n;     // number rows/columns for reduced matrix (B)
  for(r = 0; r < A->n; r++) if(drop[r] == 1) pr--;

  B = cs_spalloc(pr, pr, (pr * pr), true, false);

  cnt = 0;    // initialize running count of matrix non-zeroes for new matrix B
  cd = 0;     // initialize column decrement (records number of dropped columns)
  for(c = 0; c < A->n; c++){    // columns of A
    if(drop[c] == 1){
      cd++;
      continue;
    }

    B->p[c - cd] = cnt;

    rd = 0;   // initialize row decrement (records number of dropped rows)
    for(r = Ap[c]; r < Ap[c+1]; r++){    // rows within columns (of A)
      if(drop[Ai[r]] == 0){
        B->i[cnt] = Ai[r] - rd;
        B->x[cnt] = Ax[r];
        cnt++;
      } else rd++;

    }  // end for r (rows within columns)

  }  // end for c (columns)

  B->p[pr] = cnt;

 return(cs_done(B, NULL, NULL, 1)) ;	/* success; free workspace, return B */
}



