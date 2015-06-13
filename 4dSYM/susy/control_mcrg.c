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
  int prompt, dir, ismear, j, bl, blmax;
  int stout_step = 1;    // This might be worth reading in at some point
  double dssplaq, dstplaq, dtime, plpMod = 0.0;
  double linktr[NUMLINK], linktr_ave, linktr_width;
  complex plp = cmplx(99.0, 99.0);

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
  // Tr[Udag.U] / N
  linktr_ave = d_link(linktr, &linktr_width);
  node0_printf("FLINK");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);

  // Polyakov loop measurement and plaquette measurements
  // Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  plp = ploop(&plpMod);
  d_plaquette(&dssplaq, &dstplaq);
  node0_printf("GMES %.8g %.8g 0 %.8g %.8g ",
               plp.real, plp.imag, dssplaq, dstplaq);

  // Bosonic action (printed twice by request)
  // Might as well spit out volume average of Polyakov loop modulus
  dssplaq = d_gauge_action(NODET) / (double)volume;
  node0_printf("%.8g ", dssplaq);
  node0_printf("%.8g\n", plpMod);
  node0_printf("BACTION %.8g\n", dssplaq);

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
  for (ismear = 0; ismear <= Nstout; ismear += stout_step) {
    // For consistency, smear on unfixed links, saved here
    // Use f_U -- mom are used to store blocked links...
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->linkf[dir]), &(s->f_U[dir]));
    }

    // Calculate and print unblocked Tr[Udag.U] / N and its width
    linktr_ave = d_link(linktr, &linktr_width);
    node0_printf("BFLINK %d 0", ismear);
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", linktr[dir]);
    node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);

    // Calculate and print unblocked, plaquette determinant and widths
    blocked_plaq(ismear, 0);

#ifdef CORR
    // Konishi and SUGRA operators
    blocked_ops(ismear, 0);
#endif

    // Calculate and print unblocked Wilson loops
    blocked_rsymm(ismear, 0);

#ifdef MCRG
    // Loop over blocking levels (automatically determined)
    j = 1;
    for (bl = 1; bl <= blmax; bl++) {
      j *= 2;
      node0_printf("Blocking %d gives L = %d\n", bl, nx / j);
      block_mcrg(bl);

      // Calculate and print unblocked Tr[Udag.U] / N and its width
      linktr_ave = d_link(linktr, &linktr_width);
      node0_printf("BFLINK %d %d", ismear, bl);
      for (dir = XUP; dir < NUMLINK; dir++)
        node0_printf(" %.6g", linktr[dir]);
      node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);

      // Calculate and print blocked plaquette
      blocked_plaq(ismear, bl);

#ifdef CORR
      // Calculate and print blocked Konishi and SUGRA correlators
      blocked_ops(ismear, bl);
#endif

      // Calculate and print blocked Polyakov and Wilson loops
      blocked_ploop(ismear, bl);
      blocked_rsymm(ismear, bl);
    }
#endif

    // Here is the smearing for the next round:
    // Restore unfixed links for smearing
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->f_U[dir]), &(s->linkf[dir]));
    }
    if (ismear < Nstout) {
      node0_printf("Doing stout smearing steps %d to %d with rho=%.4g...\n",
                   ismear, ismear + stout_step, rho);

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
