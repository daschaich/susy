// -----------------------------------------------------------------
// Main procedure for N=4 SYM RG blocking
// and measurements of blocked observables including scalar correlators,
// Wilson loops and discrete R symmetry observables
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  register int i;
  register site *s;
  int prompt, dir, istout, j, bl, blmax;
  int stout_step = 1;    // This might be worth reading in at some point
  double dssplaq, dstplaq, dtime;
  complex plp = cmplx(99, 99);

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  epsilon();
  setup_PtoP();
  setup_FQ();

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Maximum number of blockings determined by smallest dimension
  if (nx < nt)
    j = nx;
  else
    j = nt;

  // Works even if we can only block down to odd j > 4
  blmax = 0;
  while (j % 2 == 0 && j > 2) {    // While j is even
    j /= 2;
    blmax++;
  }

  // Check: compute initial plaquette and bosonic action
  d_plaquette(&dssplaq, &dstplaq);
  node0_printf("START %.8g %.8g %.8g ", dssplaq, dstplaq, dssplaq + dstplaq);
  dssplaq = d_gauge_action(NODET);
  node0_printf("%.8g\n", dssplaq / (double)volume);

  // Do "local" measurements to check configuration
  // Polyakov loop measurement
  plp = ploop();

  // Tr[Udag.U] / N and plaquette measurements
  d_link(0);
  d_plaquette(&dssplaq, &dstplaq);

  // Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  node0_printf("GMES %.8g %.8g 0 %.8g %.8g ",
               plp.real, plp.imag, dssplaq, dstplaq);

  // Bosonic action (printed twice by request)
  dssplaq = d_gauge_action(NODET);
  node0_printf("%.8g\n", dssplaq / (double)volume);
  node0_printf("BACTION %.8g\n", dssplaq / (double)volume);

  // Plaquette determinant
  measure_det();

  // Monitor widths of plaquette and plaquette determinant distributions
  widths();

  // Require compilation with stout smearing enabled
#ifndef STOUT
  printf("ERROR: MCRG uses stout smearing, at least for now\n");
  terminate(1);
#endif

  // Set up P matrix for Konishi and SUGRA operators
#ifdef CORR
  setup_P();
#endif

#define MIN_PLAQ
  // Smear (after measurements) in increments of five steps up to Nstout
  for (istout = 0; istout <= Nstout; istout += stout_step) {
    // For consistency, smear on unfixed links, saved here
    // Use f_U -- mom are used to store blocked links...
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->linkf[dir]), &(s->f_U[dir]));
    }

    // Calculate and print unblocked Tr[Udag.U] / N and plaquette
    // Latter also prints unblocked plaquette determinant
    // and monitors widths of distributions
    d_link(0);
    blocked_plaq(istout, 0);

#ifdef CORR
    // Konishi and SUGRA operators
    blocked_ops(istout, 0);
#endif

    // Calculate and print unblocked Wilson loops
    blocked_rsymm(istout, 0);

#ifdef MCRG
    // Loop over blocking levels (automatically determined)
    j = 1;
    for (bl = 1; bl <= blmax; bl++) {
      j *= 2;
      node0_printf("Blocking %d gives L = %d\n", bl, nx / j);
      block_mcrg(bl);

      // Calculate and print blocked Tr[Udag.U] / N and plaquette
      d_link(bl);
      blocked_plaq(istout, bl);

#ifdef CORR
      // Calculate and print blocked Konishi and SUGRA correlators
      blocked_ops(istout, bl);
#endif

      // Calculate and print blocked Polyakov and Wilson loops
      blocked_ploop(istout, bl);
      blocked_rsymm(istout, bl);
    }
#endif

    // Here is the smearing for the next round:
    // Restore unfixed links for smearing
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->f_U[dir]), &(s->linkf[dir]));
    }
    if (istout < Nstout) {
      node0_printf("Doing stout smearing steps %d to %d with rho=%.4g...\n",
                   istout, istout + stout_step, rho);

      // Check minimum plaquette in addition to averages
      node0_printf("BEFORE ");
      d_plaquette_lcl(&dssplaq, &dstplaq);    // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g\n", dssplaq, dstplaq);

      // Overwrites s->linkf, saves original values in thin_link field
      stout_smear(stout_step, rho);
      node0_printf("AFTER  ");
      d_plaquette_lcl(&dssplaq, &dstplaq);    // Prints out MIN_PLAQ
      node0_printf(" %.8g %.8g\n", dssplaq, dstplaq);
    }
  }

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
