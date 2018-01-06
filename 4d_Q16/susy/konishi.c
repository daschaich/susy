// -----------------------------------------------------------------
// Compute and print the Konishi and SUGRA operators on each timeslice
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute traces of three bilinears using two scalar field interpolating ops
// Op 0 ("P") is log of hermitian matrix P from polar decomposition
// Op 1 ("U") is traceless part of U_a Udag_a
// Bilin 0 is Tr[PP]
// Bilin 1 is (Tr[UP] + Tr[PU]) / 2
// Bilin 2 is Tr[UU]
void compute_Ba() {
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

      // U.Udag (hermitian so trace is real)
      mult_na(&(s->link[a]), &(s->link[a]), &(Ba[1][a][i]));

      // Subtract traces
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
  // traceBB[0] is Tr[PP], traceBB[1] is Tr[UU],
  // traceBB[2] is (Tr[UP] + Tr[PU]) / 2
  FORALLDIR(a) {
    FORALLDIR(b) {
      FORALLSITES(i, s) {
        traceBB[0][a][b][i] = realtrace_nn(&(Ba[0][a][i]), &(Ba[0][b][i]));
        traceBB[2][a][b][i] = realtrace_nn(&(Ba[1][a][i]), &(Ba[1][b][i]));

        traceBB[1][a][b][i] = realtrace_nn(&(Ba[0][a][i]), &(Ba[1][b][i]));
        traceBB[1][a][b][i] += realtrace_nn(&(Ba[1][a][i]), &(Ba[0][b][i]));
        traceBB[1][a][b][i] *= 0.5;
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute and print the Konishi and SUGRA operators on each timeslice
void konishi() {
  register int i;
  register site *s;
  int a, b, t, j;
  double norm, *OK[N_K], *OS[N_K];

  // Allocate and initialize Konishi and SUGRA operators on each time slice
  for (j = 0; j < N_K; j++) {
    OK[j] = malloc(sizeof(double) * nt);
    OS[j] = malloc(sizeof(double) * nt);
    for (t = 0; t < nt; t++) {
      OK[j][t] = 0.0;
      OS[j][t] = 0.0;
    }
  }

  // Compute traces of bilinears of scalar field interpolating ops
  compute_Ba();

  // Now form the zero momentum projected operators (summing across nodes)
  FORALLSITES(i, s) {
    t = s->t;
    FORALLDIR(a) {
      for (j = 0; j < N_K; j++) {
        OK[j][t] += traceBB[j][a][a][i];    // Konishi

        // Now SUGRA, averaged over 20 off-diagonal components
        FORALLDIR(b) {
          if (a == b)
            continue;
          OS[j][t] += 0.05 * traceBB[j][a][b][i];
        }
      }
    }
  }
  // Normalization removed from site loop, followed by sum over nodes
  norm = 1.0 / (Real)(nx * ny * nz);
  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++) {
      OK[j][t] *= norm;
      OS[j][t] *= norm;
    }
    for (j = 0; j < N_K; j++) {
      g_doublesum(&(OK[j][t]));
      g_doublesum(&(OS[j][t]));
    }
  }

  // Print each operator on each time slice
  // Format: TAG  t  a  op[a]
  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++)
      node0_printf("KONISHI %d %d %.8g\n", t, j, OK[j][t]);
  }
  for (t = 0; t < nt; t++) {
    for (j = 0; j < N_K; j++)
      node0_printf("SUGRA %d %d %.8g\n", t, j, OS[j][t]);
  }

  for (j = 0; j < N_K; j++) {
    free(OK[j]);
    free(OS[j]);
  }
}
// -----------------------------------------------------------------
