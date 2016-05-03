// -----------------------------------------------------------------
// Walk around specified path of fundamental links
// dir is a list of the directions in the path, with the given length
// sign is the corresponding list of which way to go in the given dir
// That is, the negative sign means take the adjoint
// Use tempmat to accumulate linkf product along path
#include "susy_includes.h"

void path(int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j;
  matrix_f *mat;
  msg_tag *mtag;

  // Initialize tempmat with first link in path
  if (sign[0] > 0) {    // Gather from site - dir[0], no adjoint
    mtag = start_gather_site(F_OFFSET(linkf[dir[0]]), sizeof(matrix_f),
                             goffset[dir[0]] + 1, EVENANDODD, gen_pt[0]);

    wait_gather(mtag);
    FORALLSITES(i, s)
      mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmat[i]));
    cleanup_gather(mtag);
  }

  if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s)
      adjoint_f(&(s->linkf[dir[0]]), &(tempmat[i]));
  }

  // Accumulate subsequent links in product in tempmat
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_nn_f then gather from site - dir[j]
      FORALLSITES(i, s)
        mult_nn_f(&(tempmat[i]), &(s->linkf[dir[j]]), &(tempmat2[i]));

      mtag = start_gather_field(tempmat2, sizeof(matrix_f),
                                goffset[dir[j]] + 1, EVENANDODD, gen_pt[0]);

      wait_gather(mtag);
      FORALLSITES(i, s)
        mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmat[i]));
      cleanup_gather(mtag);
    }

    if (sign[j] < 0) {    // Gather from site + dir[j] then mult_na_f
      mtag = start_gather_field(tempmat, sizeof(matrix_f),
                                goffset[dir[j]], EVENANDODD, gen_pt[1]);

      wait_gather(mtag);
      FORALLSITES(i, s) {
        mat = (matrix_f *)(gen_pt[1][i]);
        mult_na_f(mat, &(s->linkf[dir[j]]), &(tempmat2[i]));
      }
      FORALLSITES(i, s)   // Don't want to overwrite tempmat too soon
        mat_copy_f(&(tempmat2[i]), &(tempmat[i]));
      cleanup_gather(mtag);
    }
  }
}
// -----------------------------------------------------------------
