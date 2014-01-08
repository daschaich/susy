// -----------------------------------------------------------------
// Evaluate the Polyakov loops using general_gathers
// Compute the Polyakov loop "at" the even sites in the first two time slices
#include "generic_includes.h"

complex ploop() {
  register int i, t;
  register site *s;
  msg_tag *tag;
  complex sum  = cmplx(0.0, 0.0), plp;
  int d[2];
  d[XUP] = 0;

  // First multiply the link on every even site by the link above it
  tag = start_gather_site(F_OFFSET(linkf[TUP]), sizeof(su3_matrix_f),
                          TUP, EVEN, gen_pt[0] );
  wait_gather(tag);
  FOREVENSITES(i, s) {
    mult_su3_nn_f(&(s->linkf[TUP]), (su3_matrix_f *)gen_pt[0][i],
                  (su3_matrix_f *)&(s->tempmat1));
  }
  cleanup_gather(tag);

  for (t = 2; t < nt; t += 2) {
    d[TUP] = t;   // Distance from which to gather
    tag = start_general_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                                    d, EVEN, gen_pt[0] );
    wait_general_gather(tag);

    // Overwrite tempmat1 on the first two time slices
    // Leave the others undisturbed so we can still gather them
    FOREVENSITES(i, s) {
      if (s->t > 1)
        continue;  // Only compute on first two slices */
      mult_su3_nn_f((su3_matrix_f *)&(s->tempmat1), (su3_matrix_f *)gen_pt[0][i],
                    (su3_matrix_f *)&(s->tempmat2));
      su3mat_copy_f((su3_matrix_f *)&(s->tempmat2), (su3_matrix_f *)&(s->tempmat1));
    }
    cleanup_general_gather(tag);
  }
  FOREVENSITES(i, s) {
    if (s->t > 1)
      continue;
    plp = trace_su3_f((su3_matrix_f *)&(s->tempmat1));
    CSUM(sum, plp);
  }

  g_complexsum(&sum);
  plp.real = sum.real / ((Real)nx);
  plp.imag = sum.imag / ((Real)nx);
  return plp;
}
// -----------------------------------------------------------------
