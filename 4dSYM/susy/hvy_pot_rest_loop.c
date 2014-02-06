// -----------------------------------------------------------------
// Modified Wilson loops for fundamental links (hardwired in path,)
// Evaluate in different spatial dirs to check rotational invariance
// This version calls path2 to compute simple rectangular loops
// with all the temporal links inverted
// so gauge fixing has no effect
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void hvy_pot_rest_loop() {
  register int i;
  register site *s;
  int t_dist, x_dist[3] = {0, 0, 0}, mu, length;
  int dir[2 * (MAX_T + MAX_X)], sign[2 * (MAX_T + MAX_X)];
  int kind[2 * (MAX_T + MAX_X)];
  double wloop;
  complex c_loop;
  su3_matrix_f tmat;

  node0_printf("hvy_pot_rest_loop: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Compute and optionally check inverse matrices
  // Temporarily store in momentum matrices
  // And store the adjoint of the inverse, it transforms like the original link does
  for (mu = XUP; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      invert(&(s->linkf[mu]), &tmat);
      su3_adjoint_f(&tmat, &(s->mom[mu]));

#ifdef DEBUG_CHECK
#define INV_TOL 1e-12
      // Check inversion
      mult_su3_an_f(&(s->mom[mu]), &(s->linkf[mu]), &tmat);
      if (1 - tmat.e[0][0].real > INV_TOL
          || tmat.e[0][0].imag > INV_TOL
          || tmat.e[0][1].real > INV_TOL
          || tmat.e[0][1].imag > INV_TOL
          || tmat.e[1][0].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL
          || 1 - tmat.e[1][1].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL) {
        printf("Link inversion fails on node%d:\n", this_node);
        dumpmat_f(&tmat);
      }
      mult_su3_na_f(&(s->linkf[mu]), &(s->mom[mu]), &tmat);
      if (1 - tmat.e[0][0].real > INV_TOL
          || tmat.e[0][0].imag > INV_TOL
          || tmat.e[0][1].real > INV_TOL
          || tmat.e[0][1].imag > INV_TOL
          || tmat.e[1][0].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL
          || 1 - tmat.e[1][1].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL) {
        printf("Link inversion fails on node%d:\n", this_node);
        dumpmat_f(&tmat);
      }
#endif
   }
  }

  // Use tempmat1 to hold loop product at each site
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    // Set up rectangular path as list of dir * sign * kind
    // For now, we only consider inverted temporal links
    for (mu = XUP; mu <= ZUP; mu++) {
      for (x_dist[mu] = 0; x_dist[mu] <= MAX_X; x_dist[mu]++) {
        if (mu > 0 && x_dist[mu] == 0)
          x_dist[mu] = 1;   // Only check (0, T) "loop" once

        length = 2 * (x_dist[mu] + t_dist);
        for (i = 0; i < x_dist[mu]; i++) {
          dir[i] = mu;
          sign[i] = 1;
          kind[i] = 0;
        }
        for (i = x_dist[mu]; i < x_dist[mu] + t_dist; i++) {
          dir[i] = TUP;
          sign[i] = 1;
          kind[i] = 1;
        }
        for (i = x_dist[mu] + t_dist; i < 2 * x_dist[mu] + t_dist; i++) {
          dir[i] = mu;
          sign[i] = -1;
          kind[i] = 0;
        }
        for (i = 2 * x_dist[mu] + t_dist; i < length; i++) {
          dir[i] = TUP;
          sign[i] = -1;
          kind[i] = 1;
        }
#ifdef DEBUG_CHECK
        node0_printf("path %d %d length %d: ", x_dist[mu], t_dist, length);
        for (i = 0; i < length; i++)
          node0_printf(" %dx%d ", dir[i] * sign[i], kind[i]);
        node0_printf("\n");
#endif

        // path2 accumulates the product in tempmat1
        path2(dir, sign, kind, length);
        wloop = 0.0;
        FORALLSITES(i, s) {
          tmat = s->tempmat1;
          c_loop = trace_su3_f(&tmat);
          wloop += c_loop.real;
        }
        g_doublesum(&wloop);
        node0_printf("PREST_LOOP %d %d %d %d %.6g\n",
                      x_dist[0], x_dist[1], x_dist[2], t_dist,
                      wloop / volume);
      } // x_dist[mu]
      x_dist[mu] = 0;
    } // mu
  } // t_dist
}
// -----------------------------------------------------------------
