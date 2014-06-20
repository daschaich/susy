// -----------------------------------------------------------------
// Measure the Konishi and SUGRA scalar field correlation functions
// Use general gathers, but combine Konishi and SUGRA into single vector
// Then only need one general gather per displacement
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Both Konishi and SUGRA correlators as functions of (x, y, z, t)
void d_correlator_r() {
  register int i;
  register site *s;
  int a, b, index, x_dist, y_dist, z_dist, t_dist;
  int len = NUMLINK * NUMLINK + 1;
  int d[4] = {0, 0, 0, 0};
  Real dum = 0.0, corr;
  Real *ops = malloc(sites_on_node * len * sizeof(Real*));
  complex ctmp;
  msg_tag *tag;
  su3_matrix_f *B[NUMLINK], tmat;

  node0_printf("d_correlator_r: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // On each site, first NUMLINK * NUMLINK components of ops are SUGRA_{ab}
  // The last component is the Konishi, which we initialize to zero
  FORALLSITES(i, s) {
    index = i * len + len - 1;
    ops[index] = 0.0;
  }

  // Allocate B_a and start calculating
  for (a = 0; a < NUMLINK; a++)
    B[a] = malloc(sites_on_node * sizeof(su3_matrix_f));

  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++)
      dum += realtrace_su3_f(&(s->linkf[a]), &(s->linkf[a]));
  }       // realtrace_su3_f(A, B) returns ReTr(Adag B)
  g_doublesum(&dum);
  dum /= (Real)(NCOL * NUMLINK * volume);

  // B_a = U_a Udag_a - Sum_{a, x} ReTr(U_a Udag_a) / (5Nc vol)
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      mult_su3_na_f(&(s->linkf[a]), &(s->linkf[a]), &(B[a][i]));
      for (b = 0; b < NCOL; b++)
        B[a][i].e[b][b].real -= dum;    // Subtract volume average
    }
  }

  // Construct the operators
  FORALLSITES(i, s) {
    for (a = 0; a < NUMLINK; a++) {
      for (b = 0; b < NUMLINK; b++) {
        mult_su3_nn_f(&(B[a][i]), &(B[b][i]), &tmat);
        ctmp = trace_su3_f(&tmat);

        // Make sure Tr(B_a * B_b) really is real
        if (abs(ctmp.imag) > IMAG_TOL) {
          printf("node%d WARNING: Im(X[%d][%d][%d]) = %.4g > %.4g)\n",
                 this_node, a, b, i, ctmp.imag, IMAG_TOL);
        }

        index = i * len + a * NUMLINK + b;
        ops[index] = ctmp.real;
        if (a == b) {
          index = i * len + len - 1;
          ops[index] += ctmp.real;   // Sum, not average
        }
      }
    }
    // Subtract trace from SUGRA operator
    for (a = 0; a < NUMLINK; a++) {
      index = i * len + len - 1;
      dum = ops[index] / (Real)NUMLINK;
      index = i * len + a * NUMLINK + a;
      ops[index] -= dum;
    }
  }
  // Done with B[a]
  for (a = 0; a < NUMLINK; a++)
    free(B[a]);

  // Construct and print correlators
  // Use general gathers, at least for now
  for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
    d[XUP] = x_dist;
    for (y_dist = 0; y_dist <= MAX_X; y_dist++) {
      d[YUP] = y_dist;
      for (z_dist = 0; z_dist <= MAX_X; z_dist++) {
        d[ZUP] = z_dist;
        for (t_dist = 0; t_dist <= MAX_T; t_dist++) {
          d[TUP] = t_dist;
          tag = start_general_gather_field(ops, len * sizeof(Real),
                                           d, EVENANDODD, gen_pt[0]);
          wait_general_gather(tag);

          // Konishi
          corr = 0.0;
          FORALLSITES(i, s) {
            index = i * len + len - 1;
            corr += ops[index] * ((Real *)gen_pt[0][i])[len - 1];
          }
          g_doublesum(&corr);
          node0_printf("CORR_K %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist, corr / volume);

          // SUGRA -- average over all 25 components
          corr = 0.0;
          FORALLSITES(i, s) {
            for (a = 0; a < len - 1; a++) {
              index = i * len + a;
              corr += ops[index] * ((Real *)gen_pt[0][i])[a];
            }
          }
          g_doublesum(&corr);
          node0_printf("CORR_S %d %d %d %d %.6g\n",
                       x_dist, y_dist, z_dist, t_dist,
                       corr / (volume * NUMLINK * NUMLINK));
          cleanup_general_gather(tag);
        } // t_dist
      } // z_dist
    } // y dist
  } // x dist
  free(ops);
}
// -----------------------------------------------------------------
