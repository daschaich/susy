// -----------------------------------------------------------------
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute traces of bilinears of scalar field interpolating ops:
// 0) traceless part of the hermitian matrix returned by the polar projection
//    (use tempmat1 as temporary storage while shifting from x+a to x)
// 1) traceless part of U_a Udag_a
void compute_Ba() {
  register int i;
  register site *s;
  int a, b, j, k, index;
  complex tc;
  Real tr;
  su3_matrix_f tmat;

  FORALLSITES(i, s) {
    // Construct scalar fields
    for (a = XUP; a < NUMLINK; a++) {
      polar(&(s->linkf[a]), &tmat, &(Ba[0][a][i]));
      mult_su3_na_f(&(s->linkf[a]), &(s->linkf[a]), &(Ba[1][a][i]));

      // Subtract traces: Both are hermitian so traces are real
      for (j = 0; j < N_B; j++) {
        tc = trace_su3_f(&(Ba[j][a][i]));
        tr = tc.real / (Real)NCOL;
        for (k = 0; k < NCOL; k++)
          Ba[j][a][i].e[k][k].real -= tr;
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Im(Tr[Ba[%d][%d][%d]]) = %.4g > %.4g\n",
                 j, a, i, fabs(tc.imag), IMAG_TOL);
        }
      }
    }
  }

  // Shift polar scalar field from x+a to x to obtain gauge-invariant traceBB
  for (a = XUP; a < NUMLINK; a++)
    shiftmat(Ba[0][a], tempmat1, goffset[a] + 1);

  // Compute traces of bilinears
  // Symmetric in a <--> b but store all to simplify SUGRA computation
  // Make sure all are purely real
  // Checked that both are gauge invariant while mixed bilinear is not
  for (a = XUP; a < NUMLINK; a++) {
    for (b = a; b < NUMLINK; b++) {
      for (j = 0; j < N_B; j++) {
        for (k = j; k < N_B; k++) {   // Skip {0, 1} = {1, 0}
          index = j + k;              // !!!
          FORALLSITES(i, s) {
            mult_su3_nn_f(&(Ba[j][a][i]), &(Ba[k][b][i]), &tmat);
            tc = trace_su3_f(&tmat);
            traceBB[index][a][b][i] = tc.real;
            traceBB[index][b][a][i] = traceBB[index][a][b][i];
            if (fabs(tc.imag) > IMAG_TOL) {
              printf("WARNING: Tr(BB[%d][%d][%d]) = (%.4g, %.4g) at site %d\n",
                     index, a, b, tc.real, tc.imag, i);
            }
          }
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure and print the Konishi and SUGRA operators on each timeslice
void d_konishi() {
  register int i;
  register site *s;
  int a, b, mu, nu, t, j;
  Real tr, norm;
  double *OK[N_K], *OS[N_K];  // Konishi and SUGRA operators on each time slice

  // Allocate and initialize Konishi and SUGRA operators
  // SUGRA will be symmetric by construction, so ignore nu < mu
  for (j = 0; j < N_K; j++) {
    OK[j] = malloc(nt * sizeof(*OK[j]));
    OS[j] = malloc(nt * sizeof(*OS[j]));
    for (t = 0; t < nt; t++) {
      OK[j][t] = 0.0;
      OS[j][t] = 0.0;
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Now form the zero momentum projected operators (summing across nodes)
  // Average SUGRA over ten components with mu <= nu
  for (a = 0; a < NUMLINK; a++) {
    FORALLSITES(i, s) {
      t = s->t;
      for (j = 0; j < N_K; j++)
        OK[j][t] += traceBB[j][a][a][i];

      for (b = 0; b < NUMLINK; b++) {
        for (j = 0; j < N_K; j++) {
          // Now SUGRA with mu--nu trace subtraction
          // Symmetric and traceless by construction so ignore nu <= mu
          for (mu = 0; mu < NDIMS; mu++) {
            for (nu = mu + 1; nu < NDIMS; nu++) {
              tr = P[mu][a] * P[nu][b] + P[nu][a] * P[mu][b];
              OS[j][t] += 0.5 * tr * traceBB[j][a][b][i];
            }
          }
        }
      }
    }
  }
  // Normalization removed from site loop, followed by sum over nodes
  norm = (Real)(nx * ny * nz);
  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++) {
      OK[j][t] /= norm;
      g_doublesum(&(OK[j][t]));

      OS[j][t] /= (6.0 * norm);
      g_doublesum(&(OS[j][t]));
    }
  }

  // Print operators
  // Subtract either ensemble average or volume average (respectively)
  for (j = 0; j < N_K; j++) {
    volK[j] = OK[j][0];
    volS[j] = OS[j][0];
    for (t = 1; t < nt; t++) {
      volK[j] += OK[j][t];
      volS[j] += OS[j][t];
    }
    volK[j] /= (double)nt;
    volS[j] /= (double)nt;
  }

  for (t = 0; t < nt; t++) {
    node0_printf("KONISHI %d", t);
    for (j = 0; j < N_K; j++)
      node0_printf(" %.16g", OK[j][t] - vevK[j]);
    for (j = 0; j < N_K; j++)
      node0_printf(" %.16g", OK[j][t] - volK[j]);
    node0_printf("\n");
  }

  for (t = 0; t < nt; t++) {
    node0_printf("SUGRA %d", t);
    for (j = 0; j < N_K; j++)
      node0_printf(" %.16g", OS[j][t]);   // Assume vanishing vev
    for (j = 0; j < N_K; j++)
      node0_printf(" %.16g", OS[j][t] - volS[j]);
    node0_printf("\n");
  }

  for (j = 0; j < N_K; j++) {
    free(OK[j]);
    free(OS[j]);
  }
}
// -----------------------------------------------------------------
