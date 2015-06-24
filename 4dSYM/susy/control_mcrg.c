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
  int smear_step = 1;    // This might be worth reading in at some point
  double dssplaq, dstplaq, dtime, plpMod = 0.0;
  double linktr[NUMLINK], linktr_ave, linktr_width;
  double link_det[NUMLINK], det_ave, det_width;
  complex plp = cmplx(99.0, 99.0);

#ifndef MCRG
  node0_printf("Don't use control_mcrg unless compiling with -DMCRG!\n");
  terminate(1);
#endif

  // Require compilation with smearing enabled
#define MIN_PLAQ
#ifndef SMEAR
  printf("ERROR: MCRG uses APE smearing, at least for now\n");
  terminate(1);
#endif

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

  // Load input and run
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
  linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
  node0_printf("FLINK");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
  node0_printf("FLINK_DET");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", link_det[dir]);
  node0_printf(" %.6g %.6g\n", det_ave, det_width);

  // Polyakov loop and plaquette measurements
  // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
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

  // Set up P matrix for Konishi and SUGRA operators
#ifdef CORR
  setup_P();
#endif

  // ---------------------------------------------------------------
  // First print unsmeared unblocked observables
  // Unsmeared unblocked Tr[Udag.U] / N and its width
  linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
  node0_printf("BFLINK 0 0");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
  node0_printf("BFLINK_DET 0 0");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", link_det[dir]);
  node0_printf(" %.6g %.6g\n", det_ave, det_width);

  // Unsmeared unblocked plaquette determinant and widths
  blocked_plaq(0, 0);

#ifdef CORR
  // Unsmeared unblocked Konishi and SUGRA operators
  blocked_ops(0, 0);
#endif

  // Unsmeared unblocked Wilson loops and R symmetry transformations
  blocked_rsymm(0, 0);
  // ---------------------------------------------------------------



  // ---------------------------------------------------------------
  // Now print unblocked (but smeared) observables
  // Can probably merge in the unsmeared checks above with a little work
  // Save unsmeared unblocked links for next blocking step
  // Use f_U -- mom is used for temporary storage by the smearing routine
  FORALLSITES(i, s) {
    for (dir = XUP; dir < NUMLINK; dir++)
      su3mat_copy_f(&(s->linkf[dir]), &(s->f_U[dir]));
  }

  // Smear (after blocking) in increments of smear_step up to Nsmear
  for (ismear = 1; ismear <= Nsmear; ismear += smear_step) {
    node0_printf("Doing %d APE smearing steps ", smear_step);
    node0_printf("(total %d) with alpha=%.4g\n", ismear, alpha);

    // Check minimum plaquette in addition to averages
    node0_printf("BEFORE ");
    blocked_plaq_lcl(ismear, 0);    // Prints out MIN_PLAQ

    // Overwrite s->linkf
    blocked_APE(smear_step, alpha, YESDET, 0);
    node0_printf("AFTER  ");
    blocked_plaq_lcl(ismear, 0);    // Prints out MIN_PLAQ

    // Smeared unblocked Tr[Udag.U] / N and its width
    linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
    node0_printf("BFLINK %d 0", ismear);
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", linktr[dir]);
    node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
    node0_printf("BFLINK_DET %d 0", ismear);
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", link_det[dir]);
    node0_printf(" %.6g %.6g\n", det_ave, det_width);

    // Smeared unblocked plaquette determinant and widths
    blocked_plaq(ismear, 0);

#ifdef CORR
    // Smeared unblocked Konishi and SUGRA operators
    blocked_ops(ismear, 0);
#endif

    // Smeared unblocked Wilson loops and R symmetry transformations
    blocked_rsymm(ismear, 0);
  }

  // Restore unsmeared unblocked links from f_U
  FORALLSITES(i, s) {
    for (dir = XUP; dir < NUMLINK; dir++)
      su3mat_copy_f(&(s->f_U[dir]), &(s->linkf[dir]));
  }
  // ---------------------------------------------------------------



  // ---------------------------------------------------------------
  // Loop over blocking levels (automatically determined)
  // Can probably merge in the unblocked checks above with a little work
  j = 1;
  for (bl = 1; bl <= blmax; bl++) {
    j *= 2;
    node0_printf("Blocking %d gives L = %d\n", bl, nx / j);
    block_mcrg(bl);

    // Unsmeared blocked Tr[Udag.U] / N and its width
    linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
    node0_printf("BFLINK 0 %d", bl);
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", linktr[dir]);
    node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
    node0_printf("BFLINK_DET 0 %d", bl);
    for (dir = XUP; dir < NUMLINK; dir++)
      node0_printf(" %.6g", link_det[dir]);
    node0_printf(" %.6g %.6g\n", det_ave, det_width);

    // Unsmeared blocked plaquette
    blocked_plaq(0, bl);

#ifdef CORR
    // Unsmeared blocked Konishi and SUGRA correlators
    blocked_ops(0, bl);
#endif

    // Unsmeared blocked Polyakov and Wilson loops
    blocked_ploop(0, bl);
    blocked_rsymm(0, bl);

    // Save unsmeared blocked links for next blocking step
    // Use f_U -- mom is used for temporary storage by the blocking routine
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->linkf[dir]), &(s->f_U[dir]));
    }

    // Smear (after blocking) in increments of smear_step up to Nsmear
    for (ismear = 1; ismear <= Nsmear; ismear += smear_step) {
      node0_printf("Doing %d APE smearing steps ", smear_step);
      node0_printf("(total %d) with alpha=%.4g\n", ismear, alpha);

      // Check minimum plaquette in addition to averages
      node0_printf("BEFORE ");
      blocked_plaq_lcl(ismear, bl);    // Prints out MIN_PLAQ

      // Overwrites s->linkf
      blocked_APE(smear_step, alpha, YESDET, bl);
      node0_printf("AFTER  ");
      blocked_plaq_lcl(ismear, bl);    // Prints out MIN_PLAQ

      // Smeared blocked Tr[Udag.U] / N and its width
      linktr_ave = d_link(linktr, &linktr_width, link_det, &det_ave, &det_width);
      node0_printf("BFLINK %d %d", ismear, bl);
      for (dir = XUP; dir < NUMLINK; dir++)
        node0_printf(" %.6g", linktr[dir]);
      node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
      node0_printf("BFLINK_DET %d %d", ismear, bl);
      for (dir = XUP; dir < NUMLINK; dir++)
        node0_printf(" %.6g", link_det[dir]);
      node0_printf(" %.6g %.6g\n", det_ave, det_width);

      // Smeared blocked plaquette determinant and widths
      blocked_plaq(ismear, bl);

#ifdef CORR
      // Smeared blocked Konishi and SUGRA operators
      blocked_ops(ismear, bl);
#endif

      // Smeared blocked Wilson loops and R symmetry transformations
      blocked_rsymm(ismear, bl);
    }

    // Restore unsmeared blocked links from f_U
    FORALLSITES(i, s) {
      for (dir = XUP; dir < NUMLINK; dir++)
        su3mat_copy_f(&(s->f_U[dir]), &(s->linkf[dir]));
    }
  }
  // ---------------------------------------------------------------

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
