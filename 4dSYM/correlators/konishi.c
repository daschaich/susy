// -----------------------------------------------------------------
// Measure the Konishi and SUGRA correlation functions
// Combine Konishi and SUGRA into single structs
#include "corr_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute Konishi and SUGRA operators
// For now just consider log-polar operator
void konishi(Kops *this_ops) {
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
    for (j = 0; j < N_K; j++) {
      // Initialize
      this_ops[i].OK[j] = 0.0;
      this_ops[i].OS[j] = 0.0;

      FORALLDIR(a) {
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
  int y_start, z_start, t_start, block, meas, this_meas;
  Real one_ov_block = 1.0 / (Real)Nblock;
  Real blockNorm = 1.0 / (Real)(Nmeas * volume);
  Real std_norm = 1.0 / (Real)(Nblock * (Nblock - 1.0));
  Real ave, err, tr;
  Kcorrs *CK = malloc(Nblock * sizeof(*CK));
  Kcorrs *CS = malloc(Nblock * sizeof(*CS));

  // Loop over all displacements
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
        // Don't need negative t_dist when x, y and z are all non-positive
        if (x_dist > 0 || y_dist > 0 || z_dist > 0)
          t_start = -MAX_X;
        else
          t_start = 0;

        // Ignore any t_dist > MAX_X even if t_dist <= MAX_T
        for (t_dist = t_start; t_dist <= MAX_X; t_dist++) {
          // Compute (x, y, z, t) correlator for each block
          for (block = 0; block < Nblock; block++) {
            // Initialize this block's correlators
            for (a = 0; a < N_K; a++) {
              for (b = 0; b < N_K; b++) {
                CK[block].C[a][b] = 0.0;
                CS[block].C[a][b] = 0.0;
              }
            }

            // Accumulate within this block
            for (meas = 0; meas < Nmeas; meas++) {
              this_meas = block * Nmeas + meas;

              // Gather ops to tempops along spatial offset, using tempops2
              FORALLSITES(i, s)
                copy_ops(&(ops[this_meas][i]), &(tempops[i]));
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
              for (j = 0; j < t_dist; j++)
                shift_ops(tempops, tempops2, goffset[TUP]);
              for (j = t_dist; j < 0; j++)
                shift_ops(tempops, tempops2, goffset[TUP] + 1);

              // N_K^2 Konishi correlators
              for (a = 0; a < N_K; a++) {
                for (b = 0; b < N_K; b++) {
                  tr = 0.0;
                  FORALLSITES(i, s)
                    tr += ops[this_meas][i].OK[a] * tempops[i].OK[b];
                  g_doublesum(&tr);
                  CK[block].C[a][b] += tr;
                }
              }

              // N_K^2 SUGRA correlators
              for (a = 0; a < N_K; a++) {
                for (b = 0; b < N_K; b++) {
                  tr = 0.0;
                  FORALLSITES(i, s)
                    tr += ops[this_meas][i].OS[a] * tempops[i].OS[b];
                  g_doublesum(&tr);
                  CS[block].C[a][b] += tr;
                }
              }
            } // Nmeas


            // Let's try to keep the size of the output files under control...
            // Print out (x, y, z, t) correlator for this block
            // Normalize by blockNorm = 1.0 / (Nmeas * volume)
            for (a = 0; a < N_K; a++) {
              for (b = 0; b < N_K; b++) {
                CK[block].C[a][b] *= blockNorm;
//                node0_printf("BLOCK_K %d %d %d %d %d %d %d %.8g\n",
//                             block, x_dist, y_dist, z_dist, t_dist,
//                             a, b, CK[block].C[a][b]);
              }
            }
            // Same as above, now for SUGRA
            for (a = 0; a < N_K; a++) {
              for (b = 0; b < N_K; b++) {
                CS[block].C[a][b] *= blockNorm;
//                node0_printf("BLOCK_S %d %d %d %d %d %d %d %.8g\n",
//                             block, x_dist, y_dist, z_dist, t_dist,
//                             a, b, CK[block].C[a][b]);
              }
            }
          } // Nblock

          // Now compute and print averages and standard errors
          // Each block is already normalized above
          for (a = 0; a < N_K; a++) {
            for (b = 0; b < N_K; b++) {
              ave = CK[0].C[a][b];             // Initialize
              for (block = 1; block < Nblock; block++)
                ave += CK[block].C[a][b];
              ave *= one_ov_block;

              // Now compute variance (square of standard error)
              tr = CK[0].C[a][b] - ave;
              err = tr * tr;                      // Initialize
              for (block = 1; block < Nblock; block++) {
                tr = CK[block].C[a][b] - ave;
                err += tr * tr;
              }
              err *= std_norm;

              // Finally print
              node0_printf("CORR_K %d %d %d %d %d %d %.6g %.4g\n",
                           x_dist, y_dist, z_dist, t_dist, a, b,
                           ave, sqrt(err));
            }
          }

          // Same as above, now for SUGRA
          for (a = 0; a < N_K; a++) {
            for (b = 0; b < N_K; b++) {
              ave = CS[0].C[a][b];                // Initialize
              for (block = 1; block < Nblock; block++)
                ave += CS[block].C[a][b];
              ave *= one_ov_block;

              // Now compute variance (square of standard error)
              tr = CS[0].C[a][b] - ave;
              err = tr * tr;                      // Initialize
              for (block = 1; block < Nblock; block++) {
                tr = CS[block].C[a][b] - ave;
                err += tr * tr;
              }
              err *= std_norm;

              // Finally print
              node0_printf("CORR_S %d %d %d %d %d %d %.6g %.4g\n",
                           x_dist, y_dist, z_dist, t_dist, a, b,
                           ave, sqrt(err));
            }
          }
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist
  free(CK);
  free(CS);
}
// ----------------------------------------------------------------
