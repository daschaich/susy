// -----------------------------------------------------------------
// Evaluate the Polyakov loop after block RG blocking steps
// Use general_gathers; lattice must be divisible by 2^block in all dirs
// Use tempmat1 and tempmat2 for temporary storage
#include "susy_includes.h"

void blocked_ploop(int Nsmear, int block) {
  register int i;
  register site *s;
  int j, bl = 2, d[4] = {0, 0, 0, 0};
  complex sum = cmplx(0.0, 0.0), plp;
  msg_tag *tag;
  su3_matrix_f *mat;

  // Allow sanity check of reproducing ploop() with this routine
  if (block <= 0)
    bl = 1;

  // Set number of links to stride, bl = 2^block
  for (j = 1; j < block; j++)
    bl *= 2;

  // Copy temporal links to tempmat1
  FORALLSITES(i, s)
    su3mat_copy_f(&(s->linkf[TUP]), &(tempmat1[i]));

  // Compute the bl-strided Polyakov loop "at" ALL the sites
  // on the first bl = 2^block timeslices
  for (j = bl; j < nt; j += bl) {
    d[TUP] = j;               // Path from which to gather
    tag = start_general_gather_field(tempmat1, sizeof(su3_matrix_f),
                                     d, EVENANDODD, gen_pt[0]);
    wait_general_gather(tag);

    // Overwrite tempmat1 on the first bl time slices
    // Leave the others undisturbed so we can still gather them
    FORALLSITES(i, s) {
      if (s->t >= bl)
        continue;
      mat = (su3_matrix_f *)gen_pt[0][i];
      mult_su3_nn_f(&(tempmat1[i]), mat, &(tempmat2[i]));
      su3mat_copy_f(&(tempmat2[i]), &(tempmat1[i]));
    }
    cleanup_general_gather(tag);
  }
  FORALLSITES(i, s) {
    if (s->t >= bl)
      continue;
    plp = trace_su3_f(&(tempmat1[i]));
    CSUM(sum, plp);
  }

  // Average all the loops we just calculated
  g_complexsum(&sum);
  plp.real = sum.real / ((Real)(nx * ny * nz * bl));
  plp.imag = sum.imag / ((Real)(nx * ny * nz * bl));
  node0_printf("BPLOOP %d %d %.8g %.8g\n", Nsmear, block, plp.real, plp.imag);
}
// -----------------------------------------------------------------
