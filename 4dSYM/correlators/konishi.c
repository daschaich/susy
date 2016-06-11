// -----------------------------------------------------------------
// Measure the Konishi and SUGRA correlation functions
// Combine Konishi and SUGRA into single structs
#include "corr_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate Konishi and SUGRA operators
// For now just consider log-polar operator
void add_konishi(Kops *this_ops) {
  register int i;
  register site *s;
  int a, b, j, k;
  Real tr;
  complex tc;
  matrix tmat;

  FORALLSITES(i, s) {
    // Construct scalar fields
    FORALLDIR(a) {
      // Log of hermitian part of polar decomposition
      polar(&(s->link[a]), &(Ba[0][a][i]), &tmat);
      matrix_log(&tmat, &(Ba[0][a][i]));

      // Subtract traces for all scalar field definitions
      for (j = 0; j < N_B; j++) {
        tc = trace(&(Ba[j][a][i]));
        tr = one_ov_N * tc.real;
        for (k = 0; k < NCOL; k++)
          Ba[j][a][i].e[k][k].real -= tr;
#ifdef DEBUG_CHECK
        // Check reality of resulting scalar fields
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr[Ba[%d][%d][%d]] = (%.4g, %.4g) is not real\n",
                 j, a, i, tc.real, tc.imag);
        }
#endif
      }
    }
  }

  // Compute traces of bilinears
  // Symmetric in a <--> b but store all to simplify SUGRA computation
  // Have checked that all are purely real and gauge invariant
  FORALLDIR(a) {
    FORALLDIR(b) {
      FORALLSITES(i, s)
        traceBB[0][a][b][i] = realtrace_nn(&(Ba[0][a][i]), &(Ba[0][b][i]));
    }
  }

  // Construct the Konishi and SUGRA operators for all scalar field definitions
  // Note: Add to this_ops passed as argument
  FORALLSITES(i, s) {
    FORALLDIR(a) {
      for (j = 0; j < N_K; j++) {
        // First Konishi
        this_ops[i].OK[j] += traceBB[j][a][a][i];

        // Now SUGRA, averaged over 20 off-diagonal components
        FORALLDIR(b) {
          if (a == b)
            continue;
          this_ops[i].OS[j] += 0.05 * traceBB[j][a][b][i];
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Simple helper function to copy operators
void copy_ops(Kops *src, Kops *dest) {
  int j;
  for (j = 0; j < N_K; j++) {
    dest->OK[j] = src->OK[j];
    dest->OS[j] = src->OS[j];
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Shift all operators without parallel transport
// The dir should come from goffset
void shift_ops(Kops *dat, Kops *temp, int dir) {
  register int i;
  register site *s;
  msg_tag *mtag;

  mtag = start_gather_field(dat, sizeof(Kops), dir, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i, s)
    copy_ops((Kops *)gen_pt[0][i], &(temp[i]));
  cleanup_gather(mtag);
  FORALLSITES(i, s)
    copy_ops(&(temp[i]), &(dat[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators as functions of r
// For gathering, all operators live in Kops structs
void correlator_r() {
  register int i;
  register site *s;
  int a, b, j, x_dist, y_dist, z_dist, t_dist;
  int y_start, z_start, t_start, this_r = -1, block, meas;
  Real one_ov_block = 1.0 / (Real)Nblock, blockNorm = 1.0 / (Real)Nmeas;
  Real std_norm = 1.0 / (Real)(Nblock * (Nblock - 1.0));
  Real ave, err, tr;
  Kcorrs **CK = malloc(Nblock * sizeof(**CK));
  Kcorrs **CS = malloc(Nblock * sizeof(**CS));

  // Set up Nblock Konishi and SUGRA correlators
  // Will be initialized within block loop below
  for (block = 0; block < Nblock; block++) {
    CK[block] = malloc(total_r * sizeof(*CK));
    CS[block] = malloc(total_r * sizeof(*CS));
  }

  // Accumulate correlators within each block
  for (block = 0; block < Nblock; block++) {
    // Initialize correlators for this block
    for (j = 0; j < total_r; j++) {
      for (a = 0; a < N_K; a++) {
        for (b = 0; b < N_K; b++) {
          CK[block][j].C[a][b] = 0.0;
          CS[block][j].C[a][b] = 0.0;
        }
      }
    }

    // Loop over measurements per block
    for (meas = 0; meas < Nmeas; meas++) {
      for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
        // Don't need negative y_dist when x_dist = 0
        if (x_dist > 0)
          y_start = -MAX_X;
        else
          y_start = 0;

        for (y_dist = y_start; y_dist <= MAX_X; y_dist++) {
          // Don't need negative z_dist when both x, y non-positive
          if (x_dist > 0 || y_dist > 0)
            z_start = -MAX_X;
          else
            z_start = 0;

          for (z_dist = z_start; z_dist <= MAX_X; z_dist++) {
            // Gather ops to tempops along spatial offset, using tempops2
            FORALLSITES(i, s)
              copy_ops(&(ops[block][i]), &(tempops[i]));
            for (j = 0; j < x_dist; j++)
              shift_ops(tempops, tempops2, goffset[XUP]);
            for (j = 0; j < y_dist; j++)
              shift_ops(tempops, tempops2, goffset[YUP]);
            for (j = y_dist; j < 0; j++)
              shift_ops(tempops, tempops2, goffset[YUP] + 1);
            for (j = 0; j < z_dist; j++)
              shift_ops(tempops, tempops2, goffset[ZUP]);
            for (j = z_dist; j < 0; j++)
              shift_ops(tempops, tempops2, goffset[ZUP] + 1);

            // Don't need negative t_dist when x, y and z are all non-positive
            // Otherwise we need to start with MAX_X shifts in the -t direction
            if (x_dist > 0 || y_dist > 0 || z_dist > 0)
              t_start = -MAX_X;
            else
              t_start = 0;
            for (j = t_start; j < 0; j++)
              shift_ops(tempops, tempops2, goffset[TUP] + 1);

            // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
            for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
              // Figure out scalar distance
              tr = A4map(x_dist, y_dist, z_dist, t_dist);
              if (tr > MAX_r - 1.0e-6) {
                // Only increment t, but still need to shift in t direction
                shift_ops(tempops, tempops2, goffset[TUP]);
                continue;
              }

              // Combine four-vectors with same scalar distance
              this_r = -1;
              for (j = 0; j < total_r; j++) {
                if (fabs(tr - lookup[j]) < 1.0e-6) {
                  this_r = j;
                  break;
                }
              }
              if (this_r < 0) {
                node0_printf("correlator_r: bad scalar distance %.4g ", tr);
                node0_printf("from displacement %d %d %d %d\n",
                             x_dist, y_dist, z_dist, t_dist);
                terminate(1);
              }

              // N_K^2 Konishi correlators
              for (a = 0; a < N_K; a++) {
                for (b = 0; b < N_K; b++) {
                  tr = 0.0;
                  FORALLSITES(i, s)
                    tr += ops[block][i].OK[a] * tempops[i].OK[b];
                  g_doublesum(&tr);
                  CK[block][this_r].C[a][b] += tr;
                }
              }

              // N_K^2 SUGRA correlators
              for (a = 0; a < N_K; a++) {
                for (b = 0; b < N_K; b++) {
                  tr = 0.0;
                  FORALLSITES(i, s)
                    tr += ops[block][i].OS[a] * tempops[i].OS[b];
                  g_doublesum(&tr);
                  CS[block][this_r].C[a][b] += tr;
                }
              }

              // As we increment t, shift in t direction
              shift_ops(tempops, tempops2, goffset[TUP]);
            } // t_dist
          } // z_dist
        } // y dist
      } // x dist
    } // Nmeas
  } // Nblock

  // Now cycle through unique scalar distances
  // Compute and print averages and standard errors
  // This is also a convenient place to normalize ops by Nmeas per block
  for (j = 0; j < total_r; j++) {
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        CK[0][this_r].C[a][b] *= blockNorm;
        ave = CK[0][this_r].C[a][b];            // Initialize
        for (block = 1; j < Nblock; block++) {
          CK[block][this_r].C[a][b] *= blockNorm;
          ave += CK[block][this_r].C[a][b];
        }
        ave *= one_ov_block;

        // Now compute variance (square of standard error)
        tr = CK[0][this_r].C[a][b] - ave;
        err = tr * tr;                          // Initialize
        for (block = 1; j < Nblock; block++) {
          tr = CK[block][this_r].C[a][b] - ave;
          err += tr * tr;
        }
        err *= std_norm;

        // Finally print
        node0_printf("CORR_K %d %.6g %d %d %.6g %.4g\n", j, lookup[j], a, b,
                     ave, sqrt(err));
      }
    }
  }

  // Same as above, now for SUGRA
  for (j = 0; j < total_r; j++) {
    for (a = 0; a < N_K; a++) {
      for (b = 0; b < N_K; b++) {
        CS[0][this_r].C[a][b] *= blockNorm;
        ave = CS[0][this_r].C[a][b];            // Initialize
        for (block = 1; j < Nblock; block++) {
          CS[block][this_r].C[a][b] *= blockNorm;
          ave += CS[block][this_r].C[a][b];
        }
        ave *= one_ov_block;

        // Now compute variance (square of standard error)
        tr = CS[0][this_r].C[a][b] - ave;
        err = tr * tr;                          // Initialize
        for (block = 1; j < Nblock; block++) {
          tr = CS[block][this_r].C[a][b] - ave;
          err += tr * tr;
        }
        err *= std_norm;

        // Finally print
        node0_printf("CORR_S %d %.6g %d %d %.6g %.4g\n", j, lookup[j], a, b,
                     ave, sqrt(err));
      }
    }
  }

  for (j = 0; j < Nblock; j++) {
    free(CK[j]);
    free(CS[j]);
  }
  free(CK);
  free(CS);
}
// ----------------------------------------------------------------
