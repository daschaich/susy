// -----------------------------------------------------------------
// Wilson loops for fundamental links (hardwired in path)
// Evaluate in different spatial dirs to check rotational invariance
// This version calls path to compute simple rectangular loops
// so gauge fixing has no effect
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void hvy_pot_loop() {
  register int i;
  register site *s;
  int t_dist, x_dist[3] = {0, 0, 0}, mu, length;
  int dir[2 * (MAX_T + MAX_X)], sign[2 * (MAX_T + MAX_X)];
  Real frac = -1.0 / (Real)NCOL;
  double wloop, detloop;
  complex det_wloop, c_loop, c1, c2, mult;
  su3_matrix_f tmat;

  node0_printf("hvy_pot_loop: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Use tempmat1 to hold loop product at each site
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    // Set up rectangular path as list of dir * sign
    for (mu = XUP; mu <= ZUP; mu++) {
      for (x_dist[mu] = 0; x_dist[mu] <= MAX_X; x_dist[mu]++) {
        if (mu > 0 && x_dist[mu] == 0)
          x_dist[mu] = 1;   // Only check (0, T) "loop" once

        length = 2 * (x_dist[mu] + t_dist);
        for (i = 0; i < x_dist[mu]; i++) {
          dir[i] = mu;
          sign[i] = 1;
        }
        for (i = x_dist[mu]; i < x_dist[mu] + t_dist; i++) {
          dir[i] = TUP;
          sign[i] = 1;
        }
        for (i = x_dist[mu] + t_dist; i < 2 * x_dist[mu] + t_dist; i++) {
          dir[i] = mu;
          sign[i] = -1;
        }
        for (i = 2 * x_dist[mu] + t_dist; i < length; i++) {
          dir[i] = TUP;
          sign[i] = -1;
        }
#ifdef DEBUG_CHECK
        node0_printf("path %d %d length %d: ", x_dist[mu], t_dist, length);
        for (i = 0; i < length; i++)
          node0_printf(" %d", dir[i] * sign[i]);
        node0_printf("\n");
#endif

        // path accumulates the product in tempmat1
        path(dir, sign, length);
        wloop = 0.0;
        detloop = 0.0;
        FORALLSITES(i, s) {
          tmat = s->tempmat1;
          c_loop = trace_su3_f(&tmat);
          wloop += c_loop.real;

          // Divide out the det raised to fractional power as
          // x^y = exp[y log x]
          det_wloop = find_det(&tmat);
          c1 = clog(&det_wloop);
          CMULREAL(c1, frac, c2);
          mult = cexp(&c2);
          CMUL(c_loop, mult, c1);
          detloop += c1.real;
        }
        g_doublesum(&wloop);
        g_doublesum(&detloop);
        node0_printf("PLOT_LOOP %d %d %d %d %.6g\n",
                     x_dist[0], x_dist[1], x_dist[2], t_dist,
                     wloop / volume);
        node0_printf("DL_LOOP   %d %d %d %d %.6g\n",
                     x_dist[0], x_dist[1], x_dist[2], t_dist,
                     detloop / volume);
      } // x_dist[mu]
      x_dist[mu] = 0;
    } // mu
  } // t_dist
}
// -----------------------------------------------------------------
