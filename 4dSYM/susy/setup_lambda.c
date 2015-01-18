// -----------------------------------------------------------------
// Set up generator matrices and epsilon^{ijklm}
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Normalization: Tr(TaTb) = -delta_ab
void setup_lambda() {
  int i, j, k, l, count;
  complex inv_sqrt = cmplx(1.0 / sqrt(2.0), 0.0);
  complex i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(2.0));

#ifdef DEBUG_CHECK
  int a;
  complex trace, tt;
  node0_printf("Computing generators for U(N)\n");
#endif

  // Make sure Lambda matrices are initialized
  for (i = 0; i < DIMF; i++)
    clear_su3mat_f(&(Lambda[i]));

  // N * (N - 1) off-diagonal SU(N) generators
  // (T^{ij, +})_{kl} = i * (de_{ki} de_{lj} + de_{kj} de_{li}) / sqrt(2)
  // (T^{ij, -})_{kl} = (de_{ki} de_{lj} - de_{kj} de_{ki}) / sqrt(2)
  // Sign in second chosen to match previous values
  count = 0;
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          if (k == i && l == j) {
            CSUM(Lambda[count].e[k][l], i_inv_sqrt);
            CSUM(Lambda[count + 1].e[k][l], inv_sqrt);
          }
          else if (k == j && l == i) {
            CSUM(Lambda[count].e[k][l], i_inv_sqrt);
            CDIF(Lambda[count + 1].e[k][l], inv_sqrt);
          }
        }
      }
      count += 2;
    }
  }
  if (count != NCOL * (NCOL - 1)) {
    node0_printf("ERROR: Wrong number of off-diagonal generators, ");
    node0_printf("%d vs. %d\n", count, NCOL * (NCOL - 1));
    terminate(1);
  }

  // N - 1 diagonal SU(N) generators
  // T^k = i * diag(1, 1, ..., -k, 0, ..., 0) / sqrt(k * (k + 1))
  for (i = 0; i < NCOL - 1; i++) {
    j = NCOL * (NCOL - 1) + i;    // Index after +/- above
    k = i + 1;
    i_inv_sqrt = cmplx(0.0, 1.0 / sqrt(k * (k + 1.0)));
    for (l = 0; l <= k; l++)
      Lambda[j].e[l][l] = i_inv_sqrt;
    CMULREAL(Lambda[j].e[k][k], -1.0 * k, Lambda[j].e[k][k]);
  }

  // U(1) generator i * I_N / sqrt(N)
  if (DIMF == NCOL * NCOL) {    // Allow SU(N) compilation for now
    i_inv_sqrt = cmplx(0, 1.0 / sqrt(NCOL));
    clear_su3mat_f(&(Lambda[DIMF - 1]));
    for (i = 0; i < NCOL; i++)
      Lambda[DIMF - 1].e[i][i] = i_inv_sqrt;
  }

#ifdef DEBUG_CHECK
  // Print Lambdas
  for (i = 0; i < DIMF; i++){
    node0_printf("Lambda[%d]\n",i);
    if (this_node == 0)
      dumpmat_f(&(Lambda[i]));
  }

  // Test group theory
  node0_printf("Check group theory ");
  node0_printf("Sum_a Lambda^a_{kl} Lambda^a_{ij} = -delta_kj delta_il\n");
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          trace = cmplx(0, 0);
          for (a = 0; a < DIMF; a++) {
            CMUL(Lambda[a].e[k][l], Lambda[a].e[i][j], tt);
            CSUM(trace, tt);
          }
          if (cabs_sq(&trace) > 1e-8)
            node0_printf("Sum_a Lambda^a_{%d%d} Lambda^a_{%d%d} = (%.4g, %.4g)\n",
                         k, j, i, l, trace.real, trace.imag);
        }
      }
    }
  }
#endif

  // Test orthogonality and compute products of Lambdas for fermion forces
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      mult_su3_nn_f(&(Lambda[i]), &(Lambda[j]), &(Lambda_prod[i][j]));
#ifdef DEBUG_CHECK
      trace = trace_su3_f(&(Lambda_prod[i][j]));
      if (trace.real * trace.real > 1e-6)
        node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n", i, j, trace.real, trace.imag);
#endif
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up five-index totally anti-symmetric tensor
// Initialize swap to avoid optimization dependence!!!
Real order(int i, int j, int k, int l, int m) {
  int seq[5] = {i, j, k, l, m};
  int swap = 1, tmp, p, permutation = 1;
  while (swap > 0) {
    swap = 0;
    for (p = 0; p < 4; p++) {
      if (seq[p] > seq[p + 1]) {
        tmp = seq[p];
        seq[p] = seq[p + 1];
        seq[p + 1] = tmp;
        swap++;
        permutation *= -1;
      }
    }
  }
  return (Real)permutation;
}


// Set up translation of (mu, nu) to linear index of anti-symmetric matrix
void setup_plaq_index() {
  int mu, nu, index;
  for (mu = 0; mu < NUMLINK; mu++) {
    plaq_index[mu][mu] = -1;
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      index = mu * (NUMLINK - 1) - mu * (mu + 1) / 2 + nu - 1;
      plaq_index[mu][nu] = index;
      plaq_index[nu][mu] = index;
    }
  }
}

void epsilon() {
  int i, j, k, l, m;
  setup_plaq_index();
  for (i = 0; i < NUMLINK; i++) {
    for (j = 0; j < NUMLINK; j++) {
      for (k = 0; k < NUMLINK; k++) {
        for (l = 0; l < NUMLINK; l++) {
          for (m = 0; m < NUMLINK; m++)
            perm[i][j][k][l][m] = 0;
        }
      }
    }
  }

  for (i = 0; i < NUMLINK; i++) {
    for (j = 0; j < NUMLINK; j++) {
      if (j == i)
        continue;
      for (k = 0; k < NUMLINK; k++) {
        if (k == j || k == i)
          continue;
        for (l = 0; l < NUMLINK; l++) {
          if (l == k || l == j || l == i)
            continue;
          for (m = 0; m < NUMLINK; m++) {
            if (m == l || m == k || m == j || m == i)
              continue;
            perm[i][j][k][l][m] = order(i, j, k, l, m);
#ifdef DEBUG_CHECK
            if (perm[i][j][k][l][m] * perm[i][j][k][l][m] > 1e-4)
              node0_printf("PERM(%d, %d, %d, %d, %d) = %.4g\n",
                           i, j, k, l, m, perm[i][j][k][l][m]);
#endif
          }
        }
      }
    }
  }
  return;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up lookup tables to simplify loops in utilities.c
// Each line is {a, b, c, d, e} for both Dbplus and Dbminus
// If counter exceeds NTERMS (amount allocated) we have a problem
void setup_PtoP() {
  int a, b, c, d, e, counter = 0;
  for (a = 0; a < NUMLINK; a++) {
    for (b = a + 1; b < NUMLINK; b++) {
      for (c = 0; c < NUMLINK; c++) {
        if (c == a || c == b)
          continue;
        for (d = 0; d < NUMLINK; d++) {
          if (d == c || d == a || d == b)
            continue;
          for (e = d + 1; e < NUMLINK; e++) {
            if (e == c || e == a || e == b)
              continue;
            DbplusPtoP_lookup[counter][0] = a;
            DbplusPtoP_lookup[counter][1] = b;
            DbplusPtoP_lookup[counter][2] = c;
            DbplusPtoP_lookup[counter][3] = d;
            DbplusPtoP_lookup[counter][4] = e;
            counter++;
          }
        }
      }
    }
  }
  if (counter > NTERMS) {
    node0_printf("ERROR: Too many terms in DbplusPtoP_lookup\n");
    terminate(1);
  }

  counter = 0;
  for (d = 0; d < NUMLINK; d++) {
    for (e = d + 1; e < NUMLINK; e++) {
      for (c = 0; c < NUMLINK; c++) {
        if (c == d || c == e)
          continue;
        for (a = 0; a < NUMLINK; a++) {
          if (a == c || a == d || a == e)
            continue;
          for (b = a + 1; b < NUMLINK; b++) {
            if (b == c || b == d || b == e)
              continue;
            DbminusPtoP_lookup[counter][0] = a;
            DbminusPtoP_lookup[counter][1] = b;
            DbminusPtoP_lookup[counter][2] = c;
            DbminusPtoP_lookup[counter][3] = d;
            DbminusPtoP_lookup[counter][4] = e;
            counter++;
          }
        }
      }
    }
  }
  if (counter > NTERMS) {
    node0_printf("ERROR: Too many terms in DbminusPtoP_lookup\n");
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up lookup tables to simplify loops in update_h.c
// Each line is {a, b, c, d, e} for Q-closed pieces in the force
// If counter exceeds NTERMS (amount allocated) we have a problem
void setup_FQ() {
  int a, b, c, d, e, counter = 0;
  for (c = 0; c < NUMLINK; c++) {
    for (d = 0; d < NUMLINK; d++) {
      if (d == c)
        continue;
      for (e = d + 1; e < NUMLINK; e++) {
        if (e == c)
          continue;
        for (a = 0; a < NUMLINK; a++) {
          if (a == d || a == e || a == c)
            continue;
          for (b = a + 1; b < NUMLINK; b++) {
            if (b == d || b == e || b == c)
              continue;
            FQ_lookup[counter][0] = a;
            FQ_lookup[counter][1] = b;
            FQ_lookup[counter][2] = c;
            FQ_lookup[counter][3] = d;
            FQ_lookup[counter][4] = e;
            counter++;
          }
        }
      }
    }
  }
  if (counter > NTERMS) {
    node0_printf("ERROR: Too many terms in FQ_lookup\n");
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up submatrix to convert from 5- to 4-component notation
//  P[NDIMS][NUMLINK] = {{t2, -t2, 0.0, 0.0, 0.0},
//                       {t6, t6, -2.0 * t6, 0.0, 0.0},
//                       {t12, t12, t12, -3.0 * t12, 0.0},
//                       {t20, t20, t20, t20, -4.0 * t20}};
// Fifth row is just all 1 / sqrt(5)
void setup_P() {
  Real t2 = 1.0 / sqrt(2),    t6 = 1.0 / sqrt(6);
  Real t12 = 1.0 / sqrt(12),  t20 = 1.0 / sqrt(20);

  P[0][0] = t2;
  P[0][1] = -t2;
  P[0][2] = 0.0;
  P[0][3] = 0.0;
  P[0][4] = 0.0;

  P[1][0] = t6;
  P[1][1] = t6;
  P[1][2] = -2.0 * t6;
  P[1][3] = 0.0;
  P[1][4] = 0.0;

  P[2][0] = t12;
  P[2][1] = t12;
  P[2][2] = t12;
  P[2][3] = -3.0 * t12;
  P[2][4] = 0.0;

  P[3][0] = t20;
  P[3][1] = t20;
  P[3][2] = t20;
  P[3][3] = t20;
  P[3][4] = -4.0 * t20;
}
// -----------------------------------------------------------------
