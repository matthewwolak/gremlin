#include "MCMCglmm.h"

cs *cs_schur(const cs *A,  int split, const cs *beta){

  int i,
      j,
      n, 
      cnt;

  n = A->n;

  cs *S11, *S22, *invS11, *S12, *S21, *muC;

  S12 = cs_spalloc (split, n-split,(n-split)*split, 1, 0);
  S11 = cs_spalloc (split, split, split*split, 1, 0);

  cnt = 0;
  for(i = split; i<n; i++){ 
    S12->p[i-split] = (i-split)*split;
    for(j = 0; j<split; j++){
      S12->i[cnt] = j;
      S12->x[cnt] = A->x[A->p[i]+j];
      cnt++;
    }
  }
  S12->p[n-split] = (n-split)*split;

  cnt = 0;
  for(i = 0; i<split; i++){ 
    S11->p[i] = i*split;
    for(j = 0; j<split; j++){
      S11->i[cnt] = j;
      S11->x[cnt] = A->x[A->p[i]+j];
      cnt++;
    }
  }
  S11->p[split] = split*split;

  invS11 = cs_inv(S11);

    S21 = cs_transpose(S12, TRUE);

    muC = cs_multiply(S21, invS11);

    cnt = 0;
    for(i = 0; i<split; i++){ 
      for(j = 0; j<(n-split); j++){
        beta->x[cnt] = muC->x[muC->p[i]+j];
        cnt++;
      }
    }

    S22 = cs_multiply(muC, S12);

    cnt = 0;

    for(i = split; i<n; i++){
      for(j = split; j<n; j++){
        S22->x[cnt] = A->x[A->p[i]+j]-S22->x[cnt];
        cnt ++;
      }
    }

    cs_spfree(S11);
    cs_spfree(invS11);
    cs_spfree(S12);
    cs_spfree(S21);
    cs_spfree(muC);

    return (cs_done (S22, NULL, NULL, 1)) ;	/* success; free workspace, return C */

}                

