// -----------------------------------------------------------------
// Evaluate the Polyakov loops using general_gathers
// Compute the Polyakov loop "at" the even sites in the first two time slices
// If PLOOPDIST dist defined (currently broken!),
// report the distribution of ploop values as a function of each x y z
#include "generic_includes.h"

complex ploop() {
  register int i, t;
  register site *s;
  msg_tag *tag;
  complex sum  = cmplx(0.0, 0.0), plp;
  int d[4];
#ifdef PLOOPDIST
  int x, y, z;
  node0_printf("PLOOPDIST currently broken!");
  terminate(1);
#endif

  d[XUP] = 0;
  d[YUP] = 0;
  d[ZUP] = 0;

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
#ifdef PLOOPDIST
    // Save result in tempmat1 for t = 0 or 1 slice
    THIS IS WRONG
    ((su3_matrix_f *)(s->tempmat1)).e[0][0] = plp;
#endif
    CSUM(sum, plp);
  }

#ifdef PLOOPDIST
  // Report ploop distribution
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      for (z = 0; z < nz; z++) {
        t = (x + y + z) % 2;
        i = node_index(x, y, z, t);
        if (node_number(x, y, z, t) != mynode())
          plp = cmplx(0.0, 0.0);
        THIS IS WRONG
        else
          plp = ((su3_matrix_f *)lattice[i].tempmat1).e[0][0];
        g_sync();
        g_complexsum(&plp);
        node0_printf("PLOOP %d %d %d %.4g %.4g\n", x, y, z, plp.real, plp.imag);
      }
    }
  }
#endif

  g_complexsum(&sum);
  plp.real = sum.real / ((Real)(nx * ny * nz));
  plp.imag = sum.imag / ((Real)(nx * ny * nz));
  return plp;
}
// -----------------------------------------------------------------
