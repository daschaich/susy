// -----------------------------------------------------------------
// Measure four-link operators related to flavor symmetry
// With unitary links, this would be a plaquette
// Our links are not unitary
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Invert single link
void invert(su3_matrix_f *in, su3_matrix_f *out) {
#if (NCOL == 1)
  complex one = cmplx(1.0, 0.0);
  CDIV(one, in->e[0][0], out->e[0][0]);
#endif

#if (NCOL == 2)
#ifndef DET
  node0_printf("invert needs compilation with -DDET\n");
  terminate(1);
#endif
  complex tc, det = find_det(in);

  // {{a, b}, {c, d}}^{-1} = {{d, -b}, {-c, a}} / det
  CDIV(in->e[1][1], det, out->e[0][0]);
  CNEGATE(in->e[0][1], tc);
  CDIV(tc, det, out->e[0][1]);
  CNEGATE(in->e[1][0], tc);
  CDIV(tc, det, out->e[1][0]);
  CDIV(in->e[0][0], det, out->e[1][1]);
#endif

#if (NCOL > 2)
  node0_printf("invert only implemented for NCOL<=2!\n");
  terminate(1);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// For now print out all {a, b}
// malloc and free temporary su3_matrix
void d_flavor() {
  register int i, dir, dir2;
  register site *s;
  double sum = 0.0;
  msg_tag *mtag0, *mtag1;
  su3_matrix_f tmat, *su3mat = malloc(sites_on_node * sizeof(*su3mat));

  if (su3mat == NULL) {
    printf("d_flavor: can't malloc su3mat\n");
    fflush(stdout);
    terminate(1);
  }

  // Compute and optionally check inverse matrices
  // Temporarily store in momentum matrices
  for (dir = XUP; dir < NUMLINK; dir++) {
    FORALLSITES(i, s) {
      invert(&(s->linkf[dir]), &(s->mom[dir]));

#ifdef DEBUG_CHECK
#define INV_TOL 1e-12
      // Check inversion
      mult_su3_nn_f(&(s->mom[dir]), &(s->linkf[dir]), &tmat);
      if (1 - tmat.e[0][0].real > INV_TOL
          || tmat.e[0][0].imag > INV_TOL
          || tmat.e[0][1].real > INV_TOL
          || tmat.e[0][1].imag > INV_TOL
          || tmat.e[1][0].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL
          || 1 - tmat.e[1][1].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL) {
        printf("Link inversion fails on node%d:\n", this_node);
        dumpmat_f(&tmat);
      }
      mult_su3_nn_f(&(s->linkf[dir]), &(s->mom[dir]), &tmat);
      if (1 - tmat.e[0][0].real > INV_TOL
          || tmat.e[0][0].imag > INV_TOL
          || tmat.e[0][1].real > INV_TOL
          || tmat.e[0][1].imag > INV_TOL
          || tmat.e[1][0].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL
          || 1 - tmat.e[1][1].real > INV_TOL
          || tmat.e[1][0].imag > INV_TOL) {
        printf("Link inversion fails on node%d:\n", this_node);
        dumpmat_f(&tmat);
      }
#endif
    }
  }

  // Construct four-link operator
  // Follow plaquette calculation, but with a couple of inverse matrices
  // Looks like the inverses make this non-symmetric under dir<-->dir2
  for (dir = XUP; dir < NUMLINK; dir++) {
    for (dir2 = XUP; dir2 < NUMLINK; dir2++) {
      if (dir2 == dir)
        continue;
      sum = 0.0;
      // gen_pt[0] is Uinv_b(x+a), gen_pt[1] is U_a(x+b)
      mtag0 = start_gather_site(F_OFFSET(mom[dir2]), sizeof(su3_matrix_f),
                                goffset[dir], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir]), sizeof(su3_matrix_f),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      // su3mat = Uinv_b(x) U_a(x)
      FORALLSITES(i, s)
        mult_su3_nn_f(&(s->mom[dir2]), &(s->linkf[dir]), &su3mat[i]);

      wait_gather(mtag0);
      wait_gather(mtag1);

      // Compute tr[Udag_a(x+b) Uinv_b(x) U_a(x) {Uinv_b(x+a)}^dag]
      FORALLSITES(i, s) {
        mult_su3_na_f(&(su3mat[i]), (su3_matrix_f *)(gen_pt[0][i]), &tmat);
        sum += (double)realtrace_su3_f((su3_matrix_f *)(gen_pt[1][i]), &tmat);
      }
      g_doublesum(&sum);

      // Average over volume and print
      sum /= (double)volume;
      node0_printf("FLAVOR %d %d %.8g\n", dir, dir2, sum);
    }
  }

  free(su3mat);
}
// -----------------------------------------------------------------
