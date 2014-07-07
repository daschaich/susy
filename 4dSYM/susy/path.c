// -----------------------------------------------------------------
// Walk around specified path of fundamental links
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dir is a list of the directions in the path, with the given length
// sign is the corresponding list of which way to go in the given dir
// That is, the negative sign means take the adjoint
// Uses tempmat1 to accumulate linkf product along path
void path(int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j;
  msg_tag *mtag;

  // Initialize tempmat1 with first link in path
  if (sign[0] > 0) {    // Gather from site - dir[0], no adjoint
    mtag = start_gather_site(F_OFFSET(linkf[dir[0]]), sizeof(su3_matrix_f),
                              goffset[dir[0]] + 1, EVENANDODD, gen_pt[0]);

    wait_gather(mtag);
    FORALLSITES(i, s)
      su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(s->tempmat1));

    cleanup_gather(mtag);
  }

  if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s)
      su3_adjoint_f(&(s->linkf[dir[0]]), &(s->tempmat1));
  }

  // Accumulate subsequent links in product in tempmat1
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_su3_nn_f then gather from site - dir[j]
      FORALLSITES(i, s)
        mult_su3_nn_f(&(s->tempmat1), &(s->linkf[dir[j]]), &(s->tempmat2));

      mtag = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix_f),
                                goffset[dir[j]] + 1, EVENANDODD, gen_pt[0]);

      wait_gather(mtag);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(s->tempmat1));

      cleanup_gather(mtag);
    }

    if (sign[j] < 0) {    // Gather from site + dir[j] then mult_su3_na_f
      mtag = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                                goffset[dir[j]], EVENANDODD, gen_pt[1]);

      wait_gather(mtag);
      FORALLSITES(i, s) {
        mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]), &(s->linkf[dir[j]]),
                      &(s->tempmat2));
      }
      FORALLSITES(i, s)   // Don't want to overwrite tempmat1 too soon
        su3mat_copy_f(&(s->tempmat2), &(s->tempmat1));

      cleanup_gather(mtag);
    }
  }
}
// -----------------------------------------------------------------
