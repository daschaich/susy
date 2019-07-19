// -----------------------------------------------------------------
// Set up 3-component epsilon tensor and gamma matrices
// Check Clifford algebra if DEBUG_CHECK #defined
#include "susy_includes.h"
//#define DEBUG_CHECK
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute tensor products of 2x2 X 2x2 and 4x4 X 2x2 matrices
// In both cases, A[i][x] * B[j][y] ends up in C[2 * i + j][2 * x + y]
void direct_prod22(int A[2][2], int B[2][2], int C[4][4]) {
  int i, j, K, x, y, Z;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      K = 2 * i + j;
      for (x = 0; x < 2; x++) {
        for (y = 0; y < 2; y++) {
          Z = 2 * x + y;
          C[K][Z] = A[i][x] * B[j][y];
        }
      }
    }
  }
}

void direct_prod42(int A[4][4], int B[2][2], int C[8][8]) {
  int i, j, K, x, y, Z;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 2; j++) {
      K = 2 * i + j;
      for (x = 0; x < 4; x++) {
        for (y = 0; y < 2; y++) {
          Z = 2 * x + y;
          C[K][Z] = A[i][x] * B[j][y];
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute inner product A[i][j] * B[j][k] * C[k][m]
// All of these are 8x8 gamma matrices
// Add/subtract to Gamma123[i][m] depending on sign
void gamma_prod(int A[8][8], int B[8][8], int C[8][8], int sign) {
  int i, j, k, m;
  for (i = 0; i < 8; i++) {
    for (m = 0; m < 8; m++) {
      for (j = 0; j < 8; j++) {
        for (k = 0; k < 8; k++)
          Gamma123[i][m] += sign * A[i][j] * B[j][k] * C[k][m];
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_gamma() {
  int i, j, k, temp44[4][4];
  int sigma1[2][2] = {{0, 1}, {1, 0}};
  int m_i_sigma2[2][2] = {{0, -1}, {1, 0}};
  int sigma3[2][2] = {{1, 0}, {0, -1}};
  int unit[2][2] = {{1, 0}, {0, 1}};

  // We'll need the 3-component epsilon tensor to set up Gamma123
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++)
        epsilon[i][j][k] = 0;
    }
  }
  epsilon[0][1][2] =  1;
  epsilon[0][2][1] = -1;
  epsilon[1][2][0] =  1;
  epsilon[2][1][0] = -1;
  epsilon[2][0][1] =  1;
  epsilon[1][0][2] = -1;

  // Direct products of Pauli matrices define the 8x8 gamma matrices
  // Gamma[0] = (-i*sigma2) X (-i*sigma2) X (-i*sigma2)
  direct_prod22(m_i_sigma2, m_i_sigma2, temp44);
  direct_prod42(temp44, m_i_sigma2, Gamma[0]);

  // Gamma[1] = sigma1 X (-i*sigma2) X unit
  direct_prod22(sigma1, m_i_sigma2, temp44);
  direct_prod42(temp44, unit, Gamma[1]);

  // Gamma[2] = sigma3 X (-i*sigma2) X unit
  direct_prod22(sigma3, m_i_sigma2, temp44);
  direct_prod42(temp44, unit, Gamma[2]);

  // Gamma[3] = (-i*sigma2) X unit x sigma1
  direct_prod22(m_i_sigma2, unit, temp44);
  direct_prod42(temp44, sigma1, Gamma[3]);

  // Gamma[4] = (-i*sigma2) X unit x sigma3
  direct_prod22(m_i_sigma2, unit, temp44);
  direct_prod42(temp44, sigma3, Gamma[4]);

  // Gamma[5] = unit x sigma1 X (-i*sigma2)
  direct_prod22(unit, sigma1, temp44);
  direct_prod42(temp44, m_i_sigma2, Gamma[5]);

  // Gamma[6] = unit x sigma3 X (-i*sigma2)
  direct_prod22(unit, sigma3, temp44);
  direct_prod42(temp44, m_i_sigma2, Gamma[6]);

  // Initialize Gamma123
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++)
      Gamma123[i][j] = 0;
  }

  // Gamma123 = epsilon[i][j][k] Gamma[i] Gamma[j] Gamma[k] / 6
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        if (epsilon[i][j][k] == 0)
          continue;

        // Add product to Gamma123
        gamma_prod(Gamma[i], Gamma[j], Gamma[k], epsilon[i][j][k]);
      }
    }
  }

  // Overly cautious to avoid integer division issues
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      if (Gamma123[i][j] == 6)
        Gamma123[i][j] = 1;
      else if (Gamma123[i][j] == -6)
        Gamma123[i][j] = -1;
      else
        Gamma123[i][j] = 0;
    }
  }

#ifdef DEBUG_CHECK
  int l, m, tg[8][8];

  // Print Gammas
  node0_printf("Computing gamma matrices\n");
  for (i = 0; i < NCHIRAL_FERMION - 1; i++){
    node0_printf("Gamma[%d]\n", i);
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      for (k = 0; k < NCHIRAL_FERMION; k++)
        node0_printf("  %d", Gamma[i][j][k]);
      node0_printf("\n");
    }
  }

  node0_printf("Gamma123\n");
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    for (k = 0; k < NCHIRAL_FERMION; k++)
      node0_printf("  %d", Gamma123[j][k]);
    node0_printf("\n");
  }

  // Test Clifford algebra
  node0_printf("Two tests of the (anti-hermitian) Clifford algebra\n");
  node0_printf("1) Gamma[i] * Gamma[i] = -I\n");
  for (i = 0; i < NCHIRAL_FERMION - 1; i++) {
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      for (k = 0; k < NCHIRAL_FERMION; k++) {
        tg[j][k] = Gamma[i][j][0] * Gamma[i][0][k];
        for (l = 1; l < NCHIRAL_FERMION; l++)
          tg[j][k] += Gamma[i][j][l] * Gamma[i][l][k];

        if (j != k && tg[j][k] != 0) {
          node0_printf("WARNING: Gamma[%d]^2_%d%d = %d\n",
                       i, j, k, tg[j][k]);
        }
      }
      if (tg[j][j] != -1) {
        node0_printf("WARNING: Gamma[%d]^2_%d%d = %d\n",
                     i, j, j, tg[j][j]);
      }
    }
  }

  node0_printf("\n2) {Gamma[i], Gamma[j]} = 0 for i != j\n");
    for (i = 0; i < NCHIRAL_FERMION - 1; i++) {
      for (j = i + 1; j < NCHIRAL_FERMION - 1; j++) {
        for (k = 0; k < NCHIRAL_FERMION; k++) {
          for (l = 0; l < NCHIRAL_FERMION; l++) {
            tg[k][l] = Gamma[i][k][0] * Gamma[j][0][l]
                       + Gamma[j][k][0] * Gamma[i][0][l];
            for (m = 1; m < NCHIRAL_FERMION; m++) {
              tg[k][l] += Gamma[i][k][m] * Gamma[j][m][l];
              tg[k][l] += Gamma[j][k][m] * Gamma[i][m][l];
            }
          if (tg[k][l] != 0) {
            node0_printf("WARNING: {Gamma[%d], Gamma[%d]}_%d%d = %d\n",
                         i, j, k, l, tg[k][l]);
          }
        }
      }
    }
  }
#endif
}
// -----------------------------------------------------------------
