// -----------------------------------------------------------------
// Update the momenta with the determinant force
// Use tempmat and staple for temporary storage
#include "susy_includes.h"

double det_force(Real eps) {
  register int i, dir, dir2;
  register site *s;
  double returnit = 0;
  complex staple_det, link_det, prod_det, tforce;
  complex *force = malloc(sites_on_node * sizeof(*force));
  matrix tmat, dlink, *mat0, *mat2;
  msg_tag *tag0 = NULL, *tag1 = NULL, *tag2 = NULL;

  // Loop over directions, update mom[dir]
  FORALLDIR(dir) {
    FORALLSITES(i, s)
      force[i] = cmplx(0.0, 0.0);

    // Loop over other directions,
    // computing force from plaquettes in the dir, dir2 plane
    FORALLDIR(dir2) {
      if (dir2 != dir) {
        // Get link[dir2] from direction dir
        tag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                 goffset[dir], EVENANDODD, gen_pt[0]);

        // Start gather for the upper staple
        tag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                 goffset[dir2], EVENANDODD, gen_pt[2]);

        // Begin the computation at the dir2DOWN point
        wait_gather(tag0);
        FORALLSITES(i, s) {
          mult_an(&(s->link[dir2]), &(s->link[dir]), &tmat);
          mult_nn(&tmat, (matrix *)gen_pt[0][i], &(tempmat[i]));
        }
        // Gather this intermediate result up to home site
        tag1 = start_gather_field(tempmat, sizeof(matrix),
                                  goffset[dir2] + 1, EVENANDODD, gen_pt[1]);

        // Begin the computation of the upper staple
        // One of the links has already been gathered
        // to compute the lower staple of the site above us in dir2
        // The plaquette is staple*U^dag due to the orientation of the gathers
        wait_gather(tag2);
        FORALLSITES(i, s) {
          mat0 = (matrix *)gen_pt[0][i];
          mat2 = (matrix *)gen_pt[2][i];
          mult_nn(&(s->link[dir2]), mat2, &tmat);
          mult_na(&tmat, mat0, &(staple[i]));

          // Now we have the upper staple -- compute its force
          // S = (det[staple U^dag] - 1) * (det[staple^dag U] - 1)
          // --> F = (det[staple U^dag] - 1) * det[staple]^* * d(det U)/dU
          //       = prod_det * staple_det^* * dlink
          staple_det = find_det(&(staple[i]));
          link_det = find_det(&(s->link[dir]));

          // prod_det = kappa_u1 * (staple_det * link_det^* - 1)
          CMUL_J(staple_det, link_det, prod_det);
          CSUM(prod_det, minus1);
          CMULREAL(prod_det, kappa_u1, prod_det);

          // force = (prod_det * staple_det^*) * dlink
          CMUL_J(prod_det, staple_det, tforce);
          CSUM(force[i], tforce);
        }

        // We have gathered up the lower staple -- compute its force
        wait_gather(tag1);
        FORALLSITES(i,s) {
          staple_det = find_det((matrix *)gen_pt[1][i]);
          link_det = find_det(&(s->link[dir]));

          // prod_det = kappa_u1 * (staple_det * link_det^* - 1)
          CMUL_J(staple_det, link_det, prod_det);
          CSUM(prod_det, minus1);
          CMULREAL(prod_det, kappa_u1, prod_det);

          // force = (prod_det * staple_det^*) * dlink
          CMUL_J(prod_det, staple_det, tforce);
          CSUM(force[i], tforce);
        }
        cleanup_gather(tag0);
        cleanup_gather(tag1);
        cleanup_gather(tag2);
      }
    }

    // Now update momenta
    FORALLSITES(i, s) {
#if (NCOL == 2 || NCOL == 3 || NCOL == 4)
      adjugate(&(s->link[dir]), &dlink);
#endif
#if (NCOL > 4)
      // Determine adjugate as determinant times inverse
      // Checked that this produces the correct results for NCOL <= 4
      invert(&(s->link[dir]), &tmat);
      link_det = find_det(&(s->link[dir]));
      c_scalar_mult_mat(&tmat, &link_det, &dlink);
#endif
      c_scalar_mult_mat(&dlink, &(force[i]), &(s->f_U[dir]));
      /* and update the momentum from the gauge force --
         dif because I computed dS/dU and the adjoint because of the way it is */
      scalar_mult_dif_adj_matrix(&(s->f_U[dir]), eps, &(s->mom[dir]));
    }
  }

  // Compute average gauge force
  FORALLSITES(i, s) {
    FORALLDIR(dir)
      returnit += realtrace(&(s->f_U[dir]), &(s->f_U[dir]));
  }
  g_doublesum(&returnit);

  free(force);
  // This will be combined with the usual gauge force terms in update_h.c
  return returnit;
}
// -----------------------------------------------------------------
