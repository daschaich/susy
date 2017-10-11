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

  // Make sure Lambda matrices are initialized
  for (i = 0; i < DIMF; i++)
    clear_mat(&(Lambda[i]));

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
    i_inv_sqrt = cmplx(0.0, sqrt(one_ov_N));
    clear_mat(&(Lambda[DIMF - 1]));
    for (i = 0; i < NCOL; i++)
      Lambda[DIMF - 1].e[i][i] = i_inv_sqrt;
  }

#ifdef DEBUG_CHECK
  int a;
  complex tc;
  matrix tmat;

  // Print Lambdas
  node0_printf("Computing generators for U(N)\n");
  for (i = 0; i < DIMF; i++){
    node0_printf("Lambda[%d]\n",i);
    if (this_node == 0)
      dumpmat(&(Lambda[i]));
  }

  // Test group theory (useful reference: arXiv:1310.5353)
  node0_printf("Check group theory ");
  node0_printf("Sum_a Lambda^a_{kl} Lambda^a_{ij} = -delta_kj delta_il\n");
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOL; k++) {
        for (l = 0; l < NCOL; l++) {
          tc.real = Lambda[0].e[k][l].real * Lambda[0].e[i][j].real
                  - Lambda[0].e[k][l].imag * Lambda[0].e[i][j].imag;
          tc.imag = Lambda[0].e[k][l].imag * Lambda[0].e[i][j].real
                  + Lambda[0].e[k][l].real * Lambda[0].e[i][j].imag;
          for (a = 1; a < DIMF; a++) {
            tc.real += Lambda[a].e[k][l].real * Lambda[a].e[i][j].real
                     - Lambda[a].e[k][l].imag * Lambda[a].e[i][j].imag;
            tc.imag += Lambda[a].e[k][l].imag * Lambda[a].e[i][j].real
                     + Lambda[a].e[k][l].real * Lambda[a].e[i][j].imag;
          }
          if (cabs_sq(&tc) > IMAG_TOL)
            node0_printf("Sum_a La^a_{%d%d} La^a_{%d%d} = (%.4g, %.4g)\n",
                         k, l, i, j, tc.real, tc.imag);
        }
      }
    }
  }

  // Test orthogonality of products of Lambdas
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      mult_nn(&(Lambda[i]), &(Lambda[j]), &tmat);
      tc = trace(&tmat);
      if (tc.real * tc.real > IMAG_TOL)
        node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n",
                     i, j, tc.real, tc.imag);
    }
  }
#endif
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
  FORALLDIR(mu) {
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
  FORALLDIR(i) {
    FORALLDIR(j) {
      FORALLDIR(k) {
        FORALLDIR(l) {
          FORALLDIR(m)
            perm[i][j][k][l][m] = 0;
        }
      }
    }
  }

  FORALLDIR(i) {
    FORALLDIR(j) {
      if (j == i)
        continue;
      FORALLDIR(k) {
        if (k == j || k == i)
          continue;
        FORALLDIR(l) {
          if (l == k || l == j || l == i)
            continue;
          FORALLDIR(m) {
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
  FORALLDIR(a) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLDIR(c) {
        if (c == a || c == b)
          continue;
        FORALLDIR(d) {
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
  FORALLDIR(d) {
    for (e = d + 1; e < NUMLINK; e++) {
      FORALLDIR(c) {
        if (c == d || c == e)
          continue;
        FORALLDIR(a) {
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
  FORALLDIR(c) {
    FORALLDIR(d) {
      if (d == c)
        continue;
      for (e = d + 1; e < NUMLINK; e++) {
        if (e == c)
          continue;
        FORALLDIR(a) {
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
