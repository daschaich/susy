// -----------------------------------------------------------------
// Evaluate the Polyakov loop after block RG blocking steps
// Use general_gathers; lattice must be divisible by 2^block in all dirs
// Use tempmat and tempmat2 for temporary storage
#include "susy_includes.h"

void blocked_ploop(int block) {
  register int i, index = node_index(0);
  register site *s;
  int j, bl = 2, d = 0;
  complex sum = cmplx(0.0, 0.0), plp;
  msg_tag *tag;

  // Allow sanity check of reproducing ploop() with this routine
  if (block <= 0)
    bl = 1;

  // Set number of links to stride, bl = 2^block
  for (j = 1; j < block; j++)
    bl *= 2;

  // Copy temporal links to tempmat
  FORALLSITES(i, s)
    mat_copy(&(s->link[TUP]), &(tempmat[i]));

  // Compute the bl-strided Polyakov loop "at" ALL the sites
  // on the first bl = 2^block timeslices
  for (j = bl; j < nt; j += bl) {
    d = j;               // Path from which to gather
    tag = start_general_gather_field(tempmat, sizeof(matrix),
                                     d, EVENANDODD, gen_pt[0]);
    wait_general_gather(tag);

    // Overwrite tempmat on the first bl time slices
    // Leave the others undisturbed so we can still gather them
    FORALLSITES(i, s) {
      if (s->t >= bl)
        continue;
      mult_nn(&(tempmat[i]), (matrix *)gen_pt[0][i], &(tempmat2[i]));
      mat_copy(&(tempmat2[i]), &(tempmat[i]));
    }
    cleanup_general_gather(tag);
  }
  FORALLSITES(i, s) {
    if (s->t >= bl)
      continue;
    plp = trace(&(tempmat[i]));
    CSUM(sum, plp);
  }

  // Average all the loops we just calculated
  g_complexsum(&sum);
  plp.real = sum.real / ((Real)bl);
  plp.imag = sum.imag / ((Real)bl);
  node0_printf("BPLOOP %d %.8g %.8g\n", block, plp.real, plp.imag);
}
// -----------------------------------------------------------------
