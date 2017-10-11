// -----------------------------------------------------------------
// For mode number computations, at all sites
// construct a simple Z2 random Twist_Fermion z_rand
#include "susy_includes.h"

void Z2source() {
  register int i, j, mu;
  register site *s;
  complex rand;

  FORALLSITES(i, s) {
    clear_TF(&(z_rand[i]));
    for (j = 0; j < DIMF; j++) {                // Site fermions
#ifdef SITERAND
      rand.real = Z2_rand_no(&(s->site_prn));
      rand.imag = Z2_rand_no(&(s->site_prn));
#else
      rand.real = Z2_rand_no(&node_prn);
      rand.imag = Z2_rand_no(&node_prn);
#endif
      c_scalar_mult_sum_mat(&(Lambda[j]), &rand, &(z_rand[i].Fsite));
      FORALLDIR(mu) {                           // Link fermions
#ifdef SITERAND
        rand.real = Z2_rand_no(&(s->site_prn));
        rand.imag = Z2_rand_no(&(s->site_prn));
#else
        rand.real = Z2_rand_no(&node_prn);
        rand.imag = Z2_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &rand, &(z_rand[i].Flink[mu]));
      }
      for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
#ifdef SITERAND
        rand.real = Z2_rand_no(&(s->site_prn));
        rand.imag = Z2_rand_no(&(s->site_prn));
#else
        rand.real = Z2_rand_no(&node_prn);
        rand.imag = Z2_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &rand, &(z_rand[i].Fplaq[mu]));
      }
    }
  }
}
// -----------------------------------------------------------------
