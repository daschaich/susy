// -----------------------------------------------------------------
// Wilson loops for polar-projected static potential using path
// Links unitarized via polar projection -- overwritten!!!
// This version computes simple rectangular loops so gauge fixing has no effect
// Use tempmat for temporary storage
#include "susy_includes.h"

void hvy_pot_polar_loop() {
  register int i;
  register site *s;
  int t_dist, x_dist, mu, length;
  int dir[2 * (MAX_T + MAX_X)], sign[2 * (MAX_T + MAX_X)];
  double polarloop;
  complex c_loop;
  matrix tmat, tmat2;

  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      // Polar projection of all links
      // To be multiplied together after projecting
      // !!! Overwrites links
      polar(&(s->link[mu]), &tmat, &tmat2);
      mat_copy(&tmat, &(s->link[mu]));
    }
  }

  node0_printf("hvy_pot_polar_loop: MAX_T = %d, MAX_X = %d\n", MAX_T, MAX_X);

  // Use tempmat to hold loop product at each site
  for (t_dist = 1; t_dist <= MAX_T; t_dist++) {
    // Set up rectangular path as list of dir * sign
    for (x_dist = 0; x_dist <= MAX_X; x_dist++) {
      length = 2 * (x_dist + t_dist);
      for (i = 0; i < x_dist; i++) {
        dir[i] = XUP;
        sign[i] = 1;
      }
      for (i = x_dist; i < x_dist + t_dist; i++) {
        dir[i] = TUP;
        sign[i] = 1;
      }
      for (i = x_dist + t_dist; i < 2 * x_dist + t_dist; i++) {
        dir[i] = XUP;
        sign[i] = -1;
      }
      for (i = 2 * x_dist + t_dist; i < length; i++) {
        dir[i] = TUP;
        sign[i] = -1;
      }
#ifdef DEBUG_CHECK
      node0_printf("path %d %d length %d: ", x_dist, t_dist, length);
      for (i = 0; i < length; i++)
        node0_printf(" %d", dir[i] * sign[i]);
      node0_printf("\n");
#endif

      // path accumulates the product in tempmat
      path(dir, sign, length);
      polarloop = 0.0;
      FORALLSITES(i, s) {
        c_loop = trace(&(tempmat[i]));
        polarloop += c_loop.real;
      }
      g_doublesum(&polarloop);
      node0_printf("PLOLAR_LOOP %d %d %.6g\n",
                   x_dist, t_dist, polarloop / volume);
    }
  }
}
// -----------------------------------------------------------------
