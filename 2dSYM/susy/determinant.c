// -----------------------------------------------------------------
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef DET
complex find_det(su3_matrix_f *Q) {
  complex det;

#if (NCOL == 1)
  det = Q->e[0][0];
#endif
#if (NCOL == 2)
  complex det2;

  CMUL(Q->e[0][0], Q->e[1][1], det);
  CMUL(Q->e[0][1], Q->e[1][0], det2);
  CSUB(det, det2, det);
#endif
#if (NCOL > 2)
  int indx[NCOL];
  Real d;
  su3_matrix_f QQ;

  su3mat_copy_f(Q, &QQ);
  ludcmp_cx(&QQ, indx, &d);
  det = cmplx(d, 0.0);
  for (i = 0; i < NCOL; i++) {
    CMUL(det, QQ.e[i][i], det2);
    det = det2;
  }
#endif
//  printf("DETER %.4g %.4g\n", det.real, det.imag);
  return det;
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef DET
void adjugate(su3_matrix_f *src, su3_matrix_f *dest) {
#if (NCOL == 1)
  dest->e[0][0] = cmplx(1.0, 0.0);
#endif
#if (NCOL == 2)
  dest->e[0][0] = src->e[1][1];
  dest->e[1][1] = src->e[0][0];
  CMULREAL(src->e[1][0], -1.0, dest->e[1][0]);
  CMULREAL(src->e[0][1], -1.0, dest->e[0][1]);
#endif
#if (NCOL > 2)
  node0_printf("Haven't coded derivative for other than 2 colors\n");
  exit(1);
#endif
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef DET
void ludcmp_cx(su3_matrix_f *a, int *indx, Real *d) {
  int i, imax, j, k;
  Real big, fdum;
  complex sum, dum, ct;
  Real vv[NCOL];

  *d = 1.0;
  for (i = 0; i < NCOL; i++) {
    big = 0.0;
    for (j = 0; j < NCOL; j++) {
      ct.real = (*a).e[i][j].real * (*a).e[i][j].real
              + (*a).e[i][j].imag * (*a).e[i][j].imag;
      if (ct.real > big)
        big = ct.real;
    }
    if (big == 0.0) {
      node0_printf("Singular matrix in routine LUDCMP");
      exit(1);
    }
    vv[i] = 1.0 / sqrt(big);
  }
  for (j = 0; j < NCOL; j++) {
    for (i = 0; i < j; i++) {
      sum = (*a).e[i][j];
      for (k = 0; k < i; k++) {
        CMUL((*a).e[i][k], (*a).e[k][j],ct);
        CSUB(sum, ct, sum);
      }
      (*a).e[i][j] = sum;
    }
    big = 0.0;
    imax = 0;
    for (i = j; i < NCOL; i++) {
      sum = (*a).e[i][j];
      for (k = 0; k < j; k++) {
        CMUL((*a).e[i][k], (*a).e[k][j],ct);
        CSUB(sum, ct, sum);
      }
      (*a).e[i][j] = sum;

      if ((fdum = vv[i] * fabs(sum.real)) >= big) {
        big = fdum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < NCOL; k++) {
        dum = (*a).e[imax][k];
        (*a).e[imax][k] = (*a).e[j][k];
        (*a).e[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (j != NCOL-1) {
      dum=(*a).e[j][j];
      for (i = j + 1; i < NCOL; i++) {
        CDIV((*a).e[i][j], dum, ct);
        (*a).e[i][j] = ct;
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef DET
void measure_det() {
  register int i, dir1, dir2;
  register site *s;
//  Real redet, imdet, replaq, implaq;
  Real norm = (Real)(NUMLINK * (NUMLINK - 1) / 2 * volume);
  complex det_link, tot_det_link, tot_det_link_sq;//, plaquette;
  su3_matrix_f tmat;
  msg_tag *mtag0, *mtag1;

  tot_det_link = cmplx(0.0, 0.0);
  tot_det_link_sq = cmplx(0.0, 0.0);
  for (dir1 = TUP; dir1 < NUMLINK; dir1++) {
    for (dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(su3_matrix_f),
                                goffset[dir1], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir1]), sizeof(su3_matrix_f),
                                goffset[dir2], EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s)
        mult_su3_an_f(&(s->linkf[dir2]), &(s->linkf[dir1]), &(s->tempmat1));

      wait_gather(mtag0);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->tempmat1), (su3_matrix_f *)(gen_pt[0][i]),
                      &(s->staple));
      }
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]), &(s->staple), &tmat);
        det_link = find_det(&tmat);
//        plaquette = trace_su3_f(&tmat);
        CSUM(tot_det_link, det_link);
        tot_det_link_sq.real += det_link.real * det_link.real;
        tot_det_link_sq.imag += det_link.imag * det_link.imag;

//        redet = det_link.real - 1.0;
//        imdet = det_link.imag;
//        replaq = plaquette.real;
//        implaq = plaquette.imag;
//        node0_printf("DELIN %.4g %.4g\n", redet * redet + imdet * imdet,
//                                          replaq * replaq + implaq * implaq);
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }

  g_complexsum(&(tot_det_link));
  g_complexsum(&(tot_det_link_sq));
  CDIVREAL(tot_det_link, norm, tot_det_link);
  CDIVREAL(tot_det_link_sq, norm, tot_det_link_sq);
  node0_printf("DET %.6g %.6g %.6g %.6g\n", tot_det_link.real,
               tot_det_link.imag, tot_det_link_sq.real, tot_det_link_sq.imag);
}
#endif
// -----------------------------------------------------------------
