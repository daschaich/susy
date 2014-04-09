// -----------------------------------------------------------------
// Evaluate the Polyakov loops using general_gathers
// Compute the Polyakov loop "at" the even sites in the first two time slices
// This version considers the (exponentiated) traceless part
#include "generic_includes.h"

complex ploop_exp() {
  register int i, t;
  register site *s;
  int d[4] = {0, 0, 0, 0};    // Path from which to gather
  Real norm = -1.0 / ((Real)NCOL);    // Will be added
  complex sum  = cmplx(0.0, 0.0), plp, tc, trace;
  msg_tag *tag;

  // First subtract the trace from each temporal link
  // Temporarily store traceless part in mom[TUP]
  FORALLSITES(i, s) {
    su3mat_copy_f(&(s->linkf[TUP]), &(s->mom[TUP]));
    trace = trace_su3_f(&(s->linkf[TUP]));
    CMULREAL(trace, norm, tc);
    c_scalar_add_diag_su3_f(&(s->mom[TUP]), &tc);
  }

  // Now multiply the link on every even site by the link above it
  tag = start_gather_site(F_OFFSET(mom[TUP]), sizeof(su3_matrix_f),
                          TUP, EVEN, gen_pt[0] );
  wait_gather(tag);
  FOREVENSITES(i, s) {
    mult_su3_nn_f(&(s->mom[TUP]), (su3_matrix_f *)gen_pt[0][i],
                  &(s->tempmat1));
  }
  cleanup_gather(tag);

  for (t = 2; t < nt; t += 2) {
    d[TUP] = t;               // Path from which to gather
    tag = start_general_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                                    d, EVEN, gen_pt[0] );
    wait_general_gather(tag);

    // Overwrite tempmat1 on the first two time slices
    // Leave the others undisturbed so we can still gather them
    FOREVENSITES(i, s) {
      if (s->t > 1)
        continue;  // Only compute on first two slices
      mult_su3_nn_f(&(s->tempmat1), (su3_matrix_f *)gen_pt[0][i],
                    &(s->tempmat2));
      su3mat_copy_f(&(s->tempmat2), &(s->tempmat1));
    }
    cleanup_general_gather(tag);
  }
  FOREVENSITES(i, s) {
    if (s->t > 1)
      continue;
    plp = trace_su3_f(&(s->tempmat1));
    CSUM(sum, plp);
  }

  g_complexsum(&sum);
  CDIVREAL(sum, ((Real)(nx * ny * nz)), tc)
  sum = cexp(&tc);
  plp.real = sum.real;
  plp.imag = sum.imag;
  return plp;
}
// -----------------------------------------------------------------
