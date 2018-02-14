// -----------------------------------------------------------------
// Set up 3-component epsilon tensor and gamma matrices
// Check Clifford algebra
#include "susy_includes.h"

void setup_gamma() {
  int i, j, k;

  // Simple 3-component epsilon tensor
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

  // Initialize gamma matrices
  for (i = 0; i < NCHIRAL_FERMION - 1; i++) {
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      for (k = 0; k < NCHIRAL_FERMION; k++)
        Gamma[i].e[j][k] = 0;
    }
  }
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    for (k = 0; k < NCHIRAL_FERMION; k++)
      Gamma123.e[j][k] = 0;
  }

  // Non-zero elements blindly copied from serial code
  Gamma[0].e[0][7] = -1;
  Gamma[0].e[1][6] =  1;
  Gamma[0].e[2][5] =  1;
  Gamma[0].e[3][4] = -1;
  Gamma[0].e[4][3] =  1;
  Gamma[0].e[5][2] = -1;
  Gamma[0].e[6][1] = -1;
  Gamma[0].e[7][0] =  1;
  Gamma[1].e[0][6] = -1;
  Gamma[1].e[1][7] = -1;
  Gamma[1].e[2][4] =  1;
  Gamma[1].e[3][5] =  1;
  Gamma[1].e[4][2] = -1;
  Gamma[1].e[5][3] = -1;
  Gamma[1].e[6][0] =  1;
  Gamma[1].e[7][1] =  1;
  Gamma[2].e[0][2] = -1;
  Gamma[2].e[1][3] = -1;
  Gamma[2].e[2][0] =  1;
  Gamma[2].e[3][1] =  1;
  Gamma[2].e[4][6] =  1;
  Gamma[2].e[5][7] =  1;
  Gamma[2].e[6][4] = -1;
  Gamma[2].e[7][5] = -1;
  Gamma[3].e[0][5] = -1;
  Gamma[3].e[1][4] = -1;
  Gamma[3].e[2][7] = -1;
  Gamma[3].e[3][6] = -1;
  Gamma[3].e[4][1] =  1;
  Gamma[3].e[5][0] =  1;
  Gamma[3].e[6][3] =  1;
  Gamma[3].e[7][2] =  1;
  Gamma[4].e[0][4] = -1;
  Gamma[4].e[1][5] =  1;
  Gamma[4].e[2][6] = -1;
  Gamma[4].e[3][7] =  1;
  Gamma[4].e[4][0] =  1;
  Gamma[4].e[5][1] = -1;
  Gamma[4].e[6][2] =  1;
  Gamma[4].e[7][3] = -1;
  Gamma[5].e[0][3] = -1;
  Gamma[5].e[1][2] =  1;
  Gamma[5].e[2][1] = -1;
  Gamma[5].e[3][0] =  1;
  Gamma[5].e[4][7] = -1;
  Gamma[5].e[5][6] =  1;
  Gamma[5].e[6][5] = -1;
  Gamma[5].e[7][4] =  1;
  Gamma[6].e[0][1] = -1;
  Gamma[6].e[1][0] =  1;
  Gamma[6].e[2][3] =  1;
  Gamma[6].e[3][2] = -1;
  Gamma[6].e[4][5] = -1;
  Gamma[6].e[5][4] =  1;
  Gamma[6].e[6][7] =  1;
  Gamma[6].e[7][6] = -1;
  Gamma123.e[0][3] =  1;
  Gamma123.e[1][2] = -1;
  Gamma123.e[2][1] = -1;
  Gamma123.e[3][0] =  1;
  Gamma123.e[4][7] =  1;
  Gamma123.e[5][6] = -1;
  Gamma123.e[6][5] = -1;
  Gamma123.e[7][4] =  1;

#ifdef DEBUG_CHECK
  int l, h;
  gamma_mat tg[NCHIRAL_FERMION - 1]; 
  
  // Print Gammas
  node0_printf("Computing gamma matrices\n");
  for (i = 0; i < NCHIRAL_FERMION - 1; i++){
    node0_printf("Gamma[%d]\n",i);
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      for (k = 0; k < NCHIRAL_FERMION; k++)
        node0_printf("  %d", Gamma[i].e[j][k]);
      node0_printf("\n");
    }
  }

  node0_printf("Gamma123\n");
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    for (k = 0; k < NCHIRAL_FERMION; k++)
      node0_printf("  %d", Gamma123.e[j][k]);
    node0_printf("\n");
  }

  // Test Clifford algebra
  node0_printf("Test Clifford algebra\n");
  node0_printf("Gamma[i] * Gamma[i] = -I\n");
  for (i = 0; i < NCHIRAL_FERMION - 1; i++) {
    for (j = 0; j < NCHIRAL_FERMION; j++) {
      for (k = 0; k < NCHIRAL_FERMION; k++) {
        tg[i].e[j][k] = 0;
        for (l = 0; l < NCHIRAL_FERMION; l++)
          tg[i].e[j][k] += Gamma[i].e[j][l] * Gamma[i].e[l][k];

        if (j != k && tg[i].e[j][k] != 0) {
          node0_printf("WARNING: Gamma[%d]^2_%d%d = %d\n",
                       i, j, k, tg[i].e[j][k]);
        }
      }
      if (tg[i].e[j][j] != -1) {
        node0_printf("WARNING: Gamma[%d]^2_%d%d = %d\n",
                     i, j, j, tg[i].e[j][j]);
      }
    }
  }

  node0_printf("Gamma[i] * Gamma[j] = 0 for i != j\n");
  for (h = 0; h < NCHIRAL_FERMION - 1; h++) {
    for (i = 0; i < NCHIRAL_FERMION - 1; i++) {
      for (j = 0; j < NCHIRAL_FERMION; j++) {
        for (k = 0; k < NCHIRAL_FERMION; k++) {
          tg[i].e[j][k] = 0;
          for (l = 0; l < NCHIRAL_FERMION; l++)
            tg[i].e[j][k] += Gamma[h].e[j][l] * Gamma[i].e[l][k];

          if (tg[i].e[j][k] != 0) {
	        node0_printf("Warning -- (Gamma[%d] * Gamma[%d])_%d%d = %d\n",
                         h, i, j, k, tg[i].e[j][k]);
          }
        }
      }
    }
  }
#endif
}
// -----------------------------------------------------------------
