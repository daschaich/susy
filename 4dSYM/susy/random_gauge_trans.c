// -----------------------------------------------------------------
// Apply random gauge transformation to the given Twist_Fermion
// in addition to the gauge links
// Run in serial, considering only a single site
// that doesn't hit a boundary
#include "susy_includes.h"

void random_gauge_trans(Twist_Fermion *TF) {
  int a, b, i, j, x = 1, y = 1, z = 1, t = 1, s = node_index(x, y, z, t);
  complex tc;
  su3_matrix_f Gmat, tmat, etamat, psimat[NUMLINK], chimat[NPLAQ];

  if (this_node != 0) {
    printf("random_gauge_trans: only implemented in serial so far\n");
    fflush(stdout);
    terminate(1);
  }
  if (nx < 4 || ny < 4 || nz < 4 || nt < 4) {
    printf("random_gauge_trans: doesn't deal with boundaries, ");
    printf("needs to be run on larger volume\n");
    fflush(stdout);
    terminate(1);
  }

  // Set up random gaussian su3_matrix_f, then unitarize it
  clear_su3mat_f(&tmat);
  for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
    tc.real = gaussian_rand_no(&(lattice[0].site_prn));
    tc.imag = gaussian_rand_no(&(lattice[0].site_prn));
#else
    tc.real = gaussian_rand_no(&(lattice[0].node_prn));
    tc.imag = gaussian_rand_no(&(lattice[0].node_prn));
#endif
    c_scalar_mult_add_su3mat_f(&tmat, &(Lambda[j]), &tc, &tmat);
  }
  polar(&tmat, &Gmat);

  // Confirm unitarity or check invariance when Gmat = I
//  mult_su3_na_f(&Gmat, &Gmat, &tmat);
//  dumpmat_f(&tmat);
//  su3mat_copy_f(&tmat, &Gmat);

  // Left side of local eta
  clear_su3mat_f(&etamat);
  // Construct G eta = sum_j eta^j G Lambda^j
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fsite.c[j];
    c_scalar_mult_add_su3mat_f(&etamat, &tmat, &tc, &etamat);
  }
  // Project out eta^j = -Tr[Lambda^j G eta]
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &etamat, &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fsite.c[j]);
  }

  // Right side of local eta
  clear_su3mat_f(&etamat);
  // Construct eta Gdag = sum_j eta^j Lambda^j Gdag
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Fsite.c[j];
    c_scalar_mult_add_su3mat_f(&etamat, &tmat, &tc, &etamat);
  }
  // Project out eta^j = -Tr[eta Gdag Lambda^j]
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&etamat, &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fsite.c[j]);
  }

  // Left side of local links and psis; right side of local chis
  for (a = XUP; a < NUMLINK; a++) {
    mult_su3_nn_f(&Gmat, &(lattice[s].linkf[a]), &tmat);
    su3mat_copy_f(&tmat, &(lattice[s].linkf[a]));

    clear_su3mat_f(&(psimat[a]));
    for (j = 0; j < DIMF; j++) {
      mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
      tc = TF[s].Flink[a].c[j];
      c_scalar_mult_add_su3mat_f(&(psimat[a]), &tmat, &tc, &(psimat[a]));
    }
    for (j = 0; j < DIMF; j++) {
      mult_su3_nn_f(&(Lambda[j]), &(psimat[a]), &tmat);
      tc = trace_su3_f(&tmat);
      CNEGATE(tc, TF[s].Flink[a].c[j]);
    }

    for (b = a + 1; b < NUMLINK; b++) {
      i = plaq_index[a][b];
      clear_su3mat_f(&(chimat[i]));
      for (j = 0; j < DIMF; j++) {
        mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
        tc = TF[s].Fplaq[i].c[j];
        c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
      }
      for (j = 0; j < DIMF; j++) {
        mult_su3_nn_f(&(chimat[i]), &(Lambda[j]), &tmat);
        tc = trace_su3_f(&tmat);
        CNEGATE(tc, TF[s].Fplaq[i].c[j]);
      }
    }
  }

  // Right side of neighboring links and psis
  // TODO: Presumably we can convert this to a loop...
  s = node_index(x - 1, y, z, t);
  mult_su3_na_f(&(lattice[s].linkf[0]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[0]));
  clear_su3mat_f(&(psimat[0]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[0].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[0]), &tmat, &tc, &(psimat[0]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[0]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Flink[0].c[j]);
  }

  s = node_index(x, y - 1, z, t);
  mult_su3_na_f(&(lattice[s].linkf[1]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[1]));
  clear_su3mat_f(&(psimat[1]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[1].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[1]), &tmat, &tc, &(psimat[1]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[1]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Flink[1].c[j]);
  }

  s = node_index(x, y, z - 1, t);
  mult_su3_na_f(&(lattice[s].linkf[2]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[2]));
  clear_su3mat_f(&(psimat[2]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[2].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[2]), &tmat, &tc, &(psimat[2]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[2]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Flink[2].c[j]);
  }

  s = node_index(x, y, z, t - 1);
  mult_su3_na_f(&(lattice[s].linkf[3]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[3]));
  clear_su3mat_f(&(psimat[3]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[3].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[3]), &tmat, &tc, &(psimat[3]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[3]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Flink[3].c[j]);
  }

  s = node_index(x + 1, y + 1, z + 1, t + 1);
  mult_su3_na_f(&(lattice[s].linkf[4]), &Gmat, &tmat);
  su3mat_copy_f(&tmat, &(lattice[s].linkf[4]));
  clear_su3mat_f(&(psimat[4]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_na_f(&(Lambda[j]), &Gmat, &tmat);
    tc = TF[s].Flink[4].c[j];
    c_scalar_mult_add_su3mat_f(&(psimat[4]), &tmat, &tc, &(psimat[4]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(psimat[4]), &(Lambda[j]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Flink[4].c[j]);
  }

  // Left side of neighboring chis
  s = node_index(x - 1, y - 1, z, t);       // 01
  i = plaq_index[0][1];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x - 1, y, z - 1, t);       // 02
  i = plaq_index[0][2];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x - 1, y, z, t - 1);       // 03
  i = plaq_index[0][3];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y + 1, z + 1, t + 1);   // 04
  i = plaq_index[0][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y - 1, z - 1, t);       // 12
  i = plaq_index[1][2];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y - 1, z, t - 1);       // 13
  i = plaq_index[1][3];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x + 1, y, z + 1, t + 1);   // 14
  i = plaq_index[1][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x, y, z - 1, t - 1);       // 23
  i = plaq_index[2][3];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x + 1, y + 1, z, t + 1);   // 24
  i = plaq_index[2][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }

  s = node_index(x + 1, y + 1, z + 1, t);   // 34
  i = plaq_index[3][4];
  clear_su3mat_f(&(chimat[i]));
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&Gmat, &(Lambda[j]), &tmat);
    tc = TF[s].Fplaq[i].c[j];
    c_scalar_mult_add_su3mat_f(&(chimat[i]), &tmat, &tc, &(chimat[i]));
  }
  for (j = 0; j < DIMF; j++) {
    mult_su3_nn_f(&(Lambda[j]), &(chimat[i]), &tmat);
    tc = trace_su3_f(&tmat);
    CNEGATE(tc, TF[s].Fplaq[i].c[j]);
  }
}
// -----------------------------------------------------------------
