// -----------------------------------------------------------------
// Apply random gauge transformation to the given Twist_Fermion
// in addition to the gauge links
// Run in serial, considering only a single site
// that doesn't hit a boundary
#include "susy_includes.h"

void random_gauge_trans(Twist_Fermion *TF) {

//edited 
      int p, q, indx;
 for (p = 0; p < NUMLINK; p++) {
    plaq_index[p][p] = -1;                                  // i,i=-1  ,  00=11=22=-1
    for (q = p + 1; q < NUMLINK; q++) {
      indx = p * (NUMLINK - 1) - p * (p + 1) / 2 + q - 1;//              01=0,02=1,12=2 
      plaq_index[p][q] = indx;
      plaq_index[q][p] = indx;                             // i,j=j,i, 01=10=0,02=20=1,12=21=2 
    }
  }

  int a, b, i, j, x = 1, t = 1, s = node_index(x, t);
  complex tc;
  matrix Gmat, tmat, etamat, psimat[NUMLINK], chimat;

  if (this_node != 0) {
    printf("random_gauge_trans: only implemented in serial so far\n");
    fflush(stdout);
    terminate(1);
  }
  if (nx < 4 || nt < 4) {
    printf("random_gauge_trans: doesn't deal with boundaries, ");
    printf("needs to be run on larger volume\n");
    fflush(stdout);
    terminate(1);
  }

  // Set up random gaussian matrix, then unitarize it
  clear_mat(&tmat);
  for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
    tc.real = gaussian_rand_no(&(lattice[0].site_prn));
    tc.imag = gaussian_rand_no(&(lattice[0].site_prn));
#else
    tc.real = gaussian_rand_no(&(lattice[0].node_prn));
    tc.imag = gaussian_rand_no(&(lattice[0].node_prn));
#endif
    c_scalar_mult_sum_mat(&(Lambda[j]), &tc, &tmat);
  }
  polar(&tmat, &Gmat);

  // Confirm unitarity or check invariance when Gmat = I
//  mult_na(&Gmat, &Gmat, &tmat);
//  dumpmat(&tmat);
//  mat_copy(&tmat, &Gmat);

  // Left side of local eta
  clear_mat(&etamat);
  // Construct G eta = sum_j eta^j G Lambda^j
  for (j = 0; j < DIMF; j++) {
    mult_nn(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fsite.c[j];
    c_scalar_mult_sum_mat(&tmat, &tc, &etamat);
  }
  // Project out eta^j = -Tr[Lambda^j G eta]
  for (j = 0; j < DIMF; j++) {
    mult_nn(&(Lambda[j]), &etamat, &tmat);
    tc = trace(&tmat);
    CNEGATE(tc, TF[s].Fsite.c[j]);
  }

  // Right side of local eta
  clear_mat(&etamat);
  // Construct eta Gdag = sum_j eta^j Lambda^j Gdag
  for (j = 0; j < DIMF; j++) {
    mult_na(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Fsite.c[j];
    c_scalar_mult_sum_mat(&tmat, &tc, &etamat);
  }
  // Project out eta^j = -Tr[eta Gdag Lambda^j]
  for (j = 0; j < DIMF; j++) {
    mult_nn(&etamat, &(Lambda[j]), &tmat);
    tc = trace(&tmat);
    CNEGATE(tc, TF[s].Fsite.c[j]);
  }

  // Left side of local links and psis; right side of local chis
  FORALLDIR(a) {
    mult_nn(&Gmat, &(lattice[s].link[a]), &tmat);
    mat_copy(&tmat, &(lattice[s].link[a]));

    clear_mat(&(psimat[a]));
    for (j = 0; j < DIMF; j++) {
      mult_nn(&Gmat, &(Lambda[j]), &tmat);
      tc = TF[s].Flink[a].c[j];
      c_scalar_mult_sum_mat(&tmat, &tc, &(psimat[a]));
    }
    for (j = 0; j < DIMF; j++) {
      mult_nn(&(Lambda[j]), &(psimat[a]), &tmat);
      tc = trace(&tmat);
      CNEGATE(tc, TF[s].Flink[a].c[j]);
    }

    for (b = a + 1; b < NUMLINK; b++) {
      clear_mat(&(chimat));                  // to be edited --- chimat[i], i=plaq_index[a][b]
      for (j = 0; j < DIMF; j++) {
        mult_na(&(Lambda[j]), &Gmat, &tmat);
        tc = TF[s].Fplaq.c[j];
        c_scalar_mult_sum_mat(&tmat, &tc, &(chimat));
      }
      for (j = 0; j < DIMF; j++) {
        mult_nn(&(chimat), &(Lambda[j]), &tmat);
        tc = trace(&tmat);
        CNEGATE(tc, TF[s].Fplaq[i].c[j]);
      }
    }
  }

  // Right side of neighboring links and psis
  // TODO: Presumably we can convert this to a loop...
  s = node_index(x - 1, t);
  mult_na(&(lattice[s].link[0]), &Gmat, &tmat);
  mat_copy(&tmat, &(lattice[s].link[0]));
  clear_mat(&(psimat[0]));
  for (j = 0; j < DIMF; j++) {
    mult_na(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[0].c[j];
    c_scalar_mult_sum_mat(&tmat, &tc, &(psimat[0]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_nn(&(psimat[0]), &(Lambda[j]), &tmat);
    tc = trace(&tmat);
    CNEGATE(tc, TF[s].Flink[0].c[j]);
  }

  s = node_index(x, t - 1);
  mult_na(&(lattice[s].link[3]), &Gmat, &tmat);
  mat_copy(&tmat, &(lattice[s].link[3]));
  clear_mat(&(psimat[3]));
  for (j = 0; j < DIMF; j++) {
    mult_na(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[3].c[j];
    c_scalar_mult_sum_mat(&tmat, &tc, &(psimat[3]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_nn(&(psimat[3]), &(Lambda[j]), &tmat);
    tc = trace(&tmat);
    CNEGATE(tc, TF[s].Flink[3].c[j]);
  }

  // Left side of neighboring chi
  s = node_index(x - 1, t - 1);
  i = plaq_index[0][3];
  clear_mat(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_nn(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_sum_mat(&tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_nn(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }
}
// -----------------------------------------------------------------
