// -----------------------------------------------------------------
// Wilson loops for static potential using path
// Print all permutations to enable tree-level improvement
// This version computes simple rectangular loops so gauge fixing has no effect
// Use tempmat for temporary storage
#include "susy_includes.h"

void hvy_pot_loop(int do_det) {
  register int i;
  register site *s;
  int t_dist, x_dist[3] = {0, 0, 0}, mu, length;
  int dir[2 * (MAX_T + MAX_X)], sign[2 * (MAX_T + MAX_X)];
  double wloop;
  complex tc;
  matrix tmat;

  node0_printf("hvy_pot_loop: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Use tempmat to hold loop product at each site
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

        // path accumulates the product in tempmat
        path(dir, sign, length);
        wloop = 0.0;
        FORALLSITES(i, s) {
          if (do_det == 1)
            det_project(&(tempmat[i]), &tmat);
          else
            mat_copy(&(tempmat[i]), &tmat);

          tc = trace(&tmat);
          wloop += tc.real;
        }
        g_doublesum(&wloop);
        if (do_det == 1) {            // Braces suppress compiler complaint
          node0_printf("DL_LOOP   ");
        }
        else
          node0_printf("PLOT_LOOP ");
        node0_printf("%d %d %d %d %.6g\n", x_dist[0], x_dist[1], x_dist[2],
                     t_dist, wloop / volume);
      } // x_dist[mu]
      x_dist[mu] = 0;
    } // mu
  } // t_dist
}
// -----------------------------------------------------------------
