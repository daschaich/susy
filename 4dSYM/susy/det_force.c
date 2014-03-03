// -----------------------------------------------------------------
// Update the momenta with the determinant force
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Cofactor of src matrix omitting row row and column col
#ifdef DET
complex cofactor(su3_matrix_f *src, int row, int col) {
  int a = 0, b = 0, i, j;
  complex submat[NCOL - 1][NCOL - 1], tc1, tc2, cof;

  // Set up submatrix omitting row row and column col
  for (i = 0; i < NCOL; i++) {
    if (i != row) {
      for (j = 0; j < NCOL; j++) {
        if (j != col) {
          set_complex_equal(&(src->e[i][j]), &(submat[a][b]));
          b++;
        }
      }
      a++;
    }
    b = 0;
  }

#if (NCOL == 3)
  // Here things are easy
  CMUL(submat[0][0], submat[1][1], tc1);
  CMUL(submat[1][0], submat[0][1], tc2);
  CSUB(tc1, tc2, cof);
#endif
#if (NCOL == 4)
  // Here things are less fun, but still not tough
  // det = 00(11*22 - 21*12) - 01(10*22 - 20*12) + 02(10*21 - 20*11)
  complex tc3, sav;
  CMUL(submat[1][1], submat[2][2], tc1);
  CMUL(submat[2][1], submat[1][2], tc2);
  CSUB(tc1, tc2, tc3);
  CMUL(submat[0][0], tc3, cof);

  CMUL(submat[1][0], submat[2][2], tc1);
  CMUL(submat[2][0], submat[1][2], tc2);
  CSUB(tc2, tc1, tc3);    // Absorb negative sign
  CMUL(submat[0][1], tc3, sav);
  CSUM(cof, sav);

  CMUL(submat[1][0], submat[2][1], tc1);
  CMUL(submat[2][0], submat[1][1], tc2);
  CSUB(tc1, tc2, tc3);
  CMUL(submat[0][2], tc3, sav);
  CSUM(cof, sav);
#endif
#if (NCOL > 4)
  node0_printf("Haven't coded cofactor for more than 4 colors\n");
  exit(1);
#endif
  return cof;
}
#endif
// -----------------------------------------------------------------





// -----------------------------------------------------------------
// Transpose of the cofactor matrix
#ifdef DET
void adjugate(su3_matrix_f *src, su3_matrix_f *dest) {
#if (NCOL == 1)
  dest->e[0][0] = cmplx(1.0, 0.0);
#endif
#if (NCOL == 2)
  dest->e[0][0] = src->e[1][1];
  dest->e[1][1] = src->e[0][0];
  CMULREAL(src->e[0][1], -1.0, dest->e[0][1]);
  CMULREAL(src->e[1][0], -1.0, dest->e[1][0]);
#endif
#if (NCOL < 5)    // Given cofactor(), this should work generically
  int i, j;

  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      dest->e[i][j] = cofactor(src, j, i);    // Note transpose!
      // Extra negative sign for 01, 03, 10, 12, 21, 23, 30 and 32
      if ((i + j) % 2 == 1)
        CNEGATE(dest->e[i][j], dest->e[i][j]);
    }
  }
#endif
#if (NCOL > 4)
  node0_printf("Haven't coded derivative for more than 3 colors\n");
  exit(1);
#endif
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef DET
double det_force(Real eps) {
  register int i, dir1, dir2;
  register site *s;
  double returnit = 0;
  complex staple_det, linkf_det, prod_det, minus1, tforce, *force;
  su3_matrix_f tmat, dlink;
  msg_tag *tag0 = NULL, *tag1 = NULL, *tag2 = NULL;

  force = malloc(sites_on_node * sizeof(*force));
  minus1 = cmplx(-1.0, 0.0);

  // Loop over directions, update mom[dir1]
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    FORALLSITES(i, s)
      force[i] = cmplx(0.0, 0.0);

    // Loop over other directions,
    // computing force from plaquettes in the dir1, dir2 plane
    for (dir2 = 0; dir2 < NUMLINK; dir2++) {
      if (dir2 != dir1) {
        // Get linkf[dir2] from direction dir1
        tag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(su3_matrix_f),
                                 goffset[dir1], EVENANDODD, gen_pt[0]);

        // Start gather for the upper staple
        tag2 = start_gather_site(F_OFFSET(linkf[dir1]), sizeof(su3_matrix_f),
                                 goffset[dir2], EVENANDODD, gen_pt[2]);

        // Begin the computation at the dir2DOWN point
        wait_gather(tag0);
        FORALLSITES(i, s) {
          mult_su3_an_f(&(s->linkf[dir2]), &(s->linkf[dir1]), &tmat);
          mult_su3_nn_f(&tmat, (su3_matrix_f *)gen_pt[0][i],
                        &(s->tempmat1));
        }
        // Gather this intermediate result up to home site
        tag1 = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
                                 goffset[dir2] + 1, EVENANDODD, gen_pt[1]);

        // Begin the computation of the upper staple
        // One of the links has already been gathered
        // to compute the lower staple of the site above us in dir2
        // The plaquette is staple*U^dag due to the orientation of the gathers
        wait_gather(tag2);
        FORALLSITES(i, s) {
          mult_su3_nn_f(&(s->linkf[dir2]), (su3_matrix_f *)gen_pt[2][i],
                        &tmat);
          mult_su3_na_f(&tmat, (su3_matrix_f *)gen_pt[0][i], &(s->staple));

          // Now we have the upper staple -- compute its force
          // S = (det[staple U^dag] - 1) * (det[staple^dag U] - 1)
          // --> F = (det[staple U^dag] -1) * det[staple]^* * d(det U)/dU
          staple_det = find_det(&(s->staple));
          linkf_det = find_det(&(s->linkf[dir1]));

          // prod_det = kappa_u1 * (staple_det * linkf_det^* - 1)
          CMUL_J(staple_det, linkf_det, prod_det);
          CSUM(prod_det, minus1);
          CMULREAL(prod_det, kappa_u1, prod_det);

          // force = (prod_det * staple_det^*) * dlink
          CMUL_J(prod_det, staple_det, tforce);
          CSUM(force[i], tforce);
        }

        // We have gathered up the lower staple -- compute its force
        wait_gather(tag1);
        FORALLSITES(i,s) {
          staple_det = find_det((su3_matrix_f *)gen_pt[1][i]);
          linkf_det = find_det(&(s->linkf[dir1]));

          // prod_det = kappa_u1 * (staple_det * linkf_det^* - 1)
          CMUL_J(staple_det, linkf_det, prod_det);
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
      adjugate(&(s->linkf[dir1]), &dlink);
      c_scalar_mult_su3mat_f(&dlink, &force[i], &(tmat));
      su3_adjoint_f(&tmat, &(s->f_U[dir1]));
      /* and update the momentum from the gauge force --
         sub because I computed dS/dU and the adjoint because of the way it is */
      scalar_mult_sub_su3_matrix_f(&(s->mom[dir1]), &(s->f_U[dir1]), eps,
                                   &(s->mom[dir1]));
    }
  }

  // Compute average gauge force
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    FORALLSITES(i, s)
      returnit += realtrace_su3_f(&(s->f_U[dir1]), &(s->f_U[dir1]));
  }
  g_doublesum(&returnit);

  free(force);
  // This will be combined with the usual gauge force terms in update_h.c
  return returnit;
}
#endif
// -----------------------------------------------------------------
