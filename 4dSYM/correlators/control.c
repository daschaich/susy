// -----------------------------------------------------------------
// Main procedure for N=4 SYM evolution and measurements
#define CONTROL
#include "corr_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i;
  register site *s;
  int prompt, j, block, meas;
  int y_dist, z_dist, t_dist;
  Real one_ov_block, blockNorm, tr;
  double ss_plaq, st_plaq, dtime, vevK[N_K], vevS[N_K];
#ifdef SMEAR
  double max_plaq = 0.0;
#endif

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Lattice volume is fixed, so we just need to set up correlators once
  // First find smallest scalar distance MAX_r allowed by MAX_X
  // (assuming MAX_T >= MAX_X)
  MAX_r = 100.0 * MAX_X;                      // To be overwritten
  for (y_dist = 0; y_dist <= MAX_X; y_dist++) {
    for (z_dist = 0; z_dist <= MAX_X; z_dist++) {
      for (t_dist = 0; t_dist <= MAX_X; t_dist++) {
        tr = A4map(MAX_X + 1, y_dist, z_dist, t_dist);
        if (tr < MAX_r)
          MAX_r = tr;
      }
    }
  }
  node0_printf("MAX = %d --> r < %.6g\n\n", MAX_X, MAX_r);

  // Now count number of unique scalar distances r < MAX_r
  // Those scalar distances themselves are recorded in the lookup array
  // The corresponding normalizations counting the number of distinct
  // four-vectors per scalar distance are recorded in the norm array
  total_r = count_points(MAX_r);

  // We can also set one_ov_block and blockNorm
  // now that we have read in Nblock and Nmeas
  one_ov_block = 1.0 / (Real)Nblock;
  blockNorm = 1.0 / (Real)Nmeas;

  // Loop over blocks, number of which is read in
  for (block = 0; block < Nblock; block++) {
    // Initialize operators for this block
    FORALLSITES(i, s) {
      for (j = 0; j < N_K; j++) {
        ops[block][i].OK[j] = 0.0;
        ops[block][i].OS[j] = 0.0;
      }
    }

    // Loop over measurements per block, also read in
    for (meas = 0; meas < Nmeas; meas++) {
      j = block * Nmeas + meas;
      startlat_p = reload_lattice(RELOAD_SERIAL, cfg[j]);

      // Check: compute initial plaquette
      plaquette(&ss_plaq, &st_plaq);
      node0_printf("START %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);

#ifdef SMEAR
      // Optionally smear
      if (smearflag == STOUT_SMEAR) {
        node0_printf("Doing %d stout smearing steps ", Nsmear);
        node0_printf("with rho=%.4g\n", alpha);
      }
      else if (smearflag == APE_SMEAR) {
        node0_printf("Doing %d polar-projected APE smearing steps ", Nsmear);
        node0_printf("with alpha=%.4g\n", alpha);
      }

      // Check minimum plaquette in addition to averages
      node0_printf("BEFORE ");
      max_plaq = local_plaquette(&ss_plaq, &st_plaq);   // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g %.8g\n", ss_plaq, st_plaq, max_plaq);

      // Overwrite s->link
      if (smearflag == STOUT_SMEAR)
        stout_smear(Nsmear, alpha);
      else if (smearflag == APE_SMEAR)
        APE_smear(Nsmear, alpha, YESDET);
      node0_printf("AFTER  ");
      max_plaq = local_plaquette(&ss_plaq, &st_plaq);      // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g %.8g\n", ss_plaq, st_plaq, max_plaq);
#endif

      // Accumulate Konishi and SUGRA operators within this block
      add_konishi(ops[block]);
      node0_printf("\n");
      fflush(stdout);
    }
  }

  // Now we can compute the vev at each site
  // Use tempops to hold vev
  // This is also a convenient place to normalize ops by Nmeas per block
  // Might as well monitor volume average of average vacuum subtractions
  for (j = 0; j < N_K; j++) {
    vevK[j] = 0.0;
    vevS[j] = 0.0;
  }
  FORALLSITES(i, s) {
    // Initialization of running sums
    for (j = 0; j < N_K; j++) {
      ops[0][i].OK[j] *= blockNorm;
      ops[0][i].OS[j] *= blockNorm;

      tempops[i].OK[j] = ops[0][i].OK[j];
      tempops[i].OS[j] = ops[0][i].OS[j];
    }
    for (block = 1; block < Nblock; block++) {
      for (j = 0; j < N_K; j++) {
        ops[block][i].OK[j] *= blockNorm;
        ops[block][i].OS[j] *= blockNorm;

        tempops[i].OK[j] += ops[block][i].OK[j];
        tempops[i].OS[j] += ops[block][i].OS[j];
      }
    }
    // Normalize vev by number of blocks
    for (j = 0; j < N_K; j++) {
      tempops[i].OK[j] *= one_ov_block;
      tempops[i].OS[j] *= one_ov_block;

      vevK[j] += tempops[i].OK[j];
      vevS[j] += tempops[i].OS[j];
    }

    // Finally subtract vev from each block
    for (block = 0; block < Nblock; block++) {
      for (j = 0; j < N_K; j++) {
        ops[block][i].OK[j] -= tempops[i].OK[j];
        ops[block][i].OS[j] -= tempops[i].OS[j];
      }
    }
  }

  // Monitor volume average of average vacuum subtractions
  for (j = 0; j < N_K; j++)
    node0_printf("VEV_K %d %.8g\n", j, vevK[j] / (Real)volume);
  for (j = 0; j < N_K; j++)
    node0_printf("VEV_S %d %.8g\n", j, vevS[j] / (Real)volume);

  // Now analyze correlators of vacuum-subtracted operators,
  // including jackknife uncertainties
  correlator_r();

  node0_printf("RUNNING COMPLETED\n");

  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
