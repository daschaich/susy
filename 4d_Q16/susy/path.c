// -----------------------------------------------------------------
// Walk around specified path of links
// dir is a list of the directions in the path, with the given length
// sign is the corresponding list of which way to go in the given dir
// That is, the negative sign means take the adjoint
// Use tempmat to accumulate link product along path
// Use tempmat2 for temporary storage
#include "susy_includes.h"

void path(int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j;
  msg_tag *mtag;

  // Initialize tempmat with first link in path
  if (sign[0] > 0) {    // Gather from site - dir[0], no adjoint
    mtag = start_gather_site(F_OFFSET(link[dir[0]]), sizeof(matrix),
                             goffset[dir[0]] + 1, EVENANDODD, gen_pt[0]);

    wait_gather(mtag);
    FORALLSITES(i, s)
      mat_copy((matrix *)(gen_pt[0][i]), &(tempmat[i]));
    cleanup_gather(mtag);
  }

  if (sign[0] < 0) {    // Take adjoint, no gather
    FORALLSITES(i, s)
      adjoint(&(s->link[dir[0]]), &(tempmat[i]));
  }

  // Accumulate subsequent links in product in tempmat
  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {    // mult_nn then gather from site - dir[j]
      FORALLSITES(i, s)
        mult_nn(&(tempmat[i]), &(s->link[dir[j]]), &(tempmat2[i]));

      mtag = start_gather_field(tempmat2, sizeof(matrix),
                                goffset[dir[j]] + 1, EVENANDODD, gen_pt[0]);

      wait_gather(mtag);
      FORALLSITES(i, s)
        mat_copy((matrix *)(gen_pt[0][i]), &(tempmat[i]));
      cleanup_gather(mtag);
    }

    if (sign[j] < 0) {    // Gather from site + dir[j] then mult_na
      mtag = start_gather_field(tempmat, sizeof(matrix),
                                goffset[dir[j]], EVENANDODD, gen_pt[0]);

      // Be careful about overwriting tempmat;
      // gen_pt may just point to it for on-node "gathers"
      wait_gather(mtag);
      FORALLSITES(i, s)
        mult_na((matrix *)(gen_pt[0][i]), &(s->link[dir[j]]), &(tempmat2[i]));
      cleanup_gather(mtag);
      FORALLSITES(i, s)
        mat_copy(&(tempmat2[i]), &(tempmat[i]));
    }
  }
}
// -----------------------------------------------------------------
