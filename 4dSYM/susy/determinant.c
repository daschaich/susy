// -----------------------------------------------------------------
#include "susy_includes.h"
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
  int i, indx[NCOL];
  Real d;
  complex det2;
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
  for (dir1 = YUP; dir1 < NUMLINK; dir1++) {
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



// -----------------------------------------------------------------
// Measure density of monopole world lines in non-diagonal cubes
#ifdef DET
void monopole() {
  register int i, dir1, dir2;
  register site *s;
  int a, b, c, d, ip2;
  int total_mono_p[4], total_mono_m[4], total, total_abs;
  int *mono[NUMLINK][NUMLINK], *charge[NUMLINK];
  Real *phase[NUMLINK], p2, p3, total_phase, permm;
  double threePI = 3.0 * PI;
  double fivePI = 5.0 * PI;
  double sevenPI = 7.0 * PI;
  complex det_link;
  msg_tag *mtag0, *mtag1;

  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    charge[dir1] = malloc(sites_on_node * sizeof(int));
    phase[dir1] = malloc(sites_on_node * sizeof(Real));
  }
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    for (dir2 = 0; dir2 < NUMLINK; dir2++)
      mono[dir1][dir2] = malloc(sites_on_node * sizeof(int));
  }

  // First extract the U(1) part of the link
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    FORALLSITES(i, s) {
      det_link = find_det(&(s->linkf[dir1]));
      phase[dir1][i] = atan2(det_link.imag, det_link.real);
//      printf("XXX (%d, %d, %d, %d)[%d] %.6g %.6g %.6g\n",
//             s->x, s->y, s->z, s->t, dir1,
//             det_link.real, det_link.imag, phase[dir1][i]);
    }
  }

  // Next find the number of strings
  // penetrating the face of every plaquette
  for (dir1 = YUP; dir1 <= TUP; dir1++) {
    for (dir2 = XUP; dir2 < dir1; dir2++) {
      mtag0 = start_gather_field(phase[dir2], sizeof(Real),
                                 goffset[dir1], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_field(phase[dir1], sizeof(Real),
                                 goffset[dir2], EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        p2 = *((Real *)(gen_pt[0][i]));
        p3 = *((Real *)(gen_pt[1][i]));
        total_phase = phase[dir1][i] + p2 - p3 - phase[dir2][i];
//        printf("YYY (%d, %d, %d, %d)[%d,%d] %.6g %.6g %.6g %.6g %.6g\n",
//               s->x, s->y, s->z, s->t, dir1, dir2,
//               phase[dir1][i], phase[dir2][i], p2, p3, total_phase);

        if (abs(total_phase) < PI)
          mono[dir1][dir2][i] = 0;
        else if (total_phase >= PI && total_phase < threePI)
          mono[dir1][dir2][i] = 1;
        else if (total_phase >= threePI && total_phase < fivePI)
          mono[dir1][dir2][i] = 2;
        else if (total_phase >= fivePI && total_phase < sevenPI)
          mono[dir1][dir2][i] = 3;
        else if (total_phase <= -PI && total_phase > -threePI)
          mono[dir1][dir2][i] = -1;
        else if (total_phase <= -threePI && total_phase > -fivePI)
          mono[dir1][dir2][i] = -2;
        else if (total_phase <= -fivePI && total_phase > -sevenPI)
          mono[dir1][dir2][i] = -3;
        else {
          printf("monopole: total_phase %.4g out of bounds on node%d\n",
                 total_phase, this_node);
          terminate(1);
        }
        mono[dir2][dir1][i] = -mono[dir1][dir2][i];

//        if (i == 0) {
//          printf("ZZZ (%d, %d, %d, %d)[%d, %d] ",
//                 s->x, s->y, s->z, s->t, dir1, dir2);
//          printf("%.6g - %.6g + %.6g - %.6g = %.6g --> %d\n",
//                 phase[dir1][i], phase[dir2][i], p2, p3,
//                 total_phase, mono[dir1][dir2][i]);
//        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }

  // We have the number of strings penetrating every plaquette
  // Now tie these together into cubes,
  // using the first four components of the 5d epsilon
  for (a = 0; a < NUMLINK - 1; a++) {
    FORALLSITES(i, s)
      charge[a][i] = 0;

    for (b = 0; b < NUMLINK - 1; b++) {
      if (b != a) {
        for (c = 0; c < NUMLINK - 1; c++) {
          if (c == a || c == b)
            continue;
          for (d = 0; d < NUMLINK - 1; d++) {
            if (d == c || d == a || d == b)
              continue;

            permm = perm[a][b][c][d][NUMLINK - 1];
//            printf("CHECK perm[%d][%d][%d][%d] = %.1g \n",
//                   a, b, c, d, permm);

            mtag0 = start_gather_field(mono[c][d], sizeof(int),
                                       goffset[b], EVENANDODD, gen_pt[0]);

            wait_gather(mtag0);
            FORALLSITES(i, s) {
              ip2 = *((int *)(gen_pt[0][i]));
              if (permm > 0.0)
                charge[a][i] += (mono[c][d][i] - ip2);
              if (permm < 0.0)
                charge[a][i] -= (mono[c][d][i] - ip2);
              //if (s->x ==0 && s->y==0 && s->z==0 && s->t == 0)
              //printf("a %d b %d c %d d %d %d %d\n",a,b,c,d,mono[c][d][i],ip2);
            }
            cleanup_gather(mtag0);
          }
        }
      }
    }
    FORALLSITES(i, s)       // Normalize for epsilon loops double-counting
      charge[a][i] /= 2.0;  // Should end up an integer
  }
//  FORALLSITES(i, s) {
//    for (dir1 = XUP; dir1 < NUMLINK - 1; dir1++) {
//      if (charge[dir1][i] != 0) {
//        printf("QQQ (%d, %d, %d, %d)[%d] %d\n",
//               s->x, s->y, s->z, s->t, dir1, charge[dir1][i]);
//      }
//    }
//  }

  // Finally accumulate and print global quantities
  total = 0;
  total_abs = 0;
  for (dir1 = XUP; dir1 < NUMLINK - 1; dir1++) {
    total_mono_p[dir1] = 0;
    total_mono_m[dir1] = 0;
  }
  FORALLSITES(i, s) {
    for (dir1 = XUP; dir1 < NUMLINK - 1; dir1++) {
      if (charge[dir1][i] > 0)
        total_mono_p[dir1] += charge[dir1][i];
      if (charge[dir1][i] < 0)
        total_mono_m[dir1] += charge[dir1][i];
    }
    // This is curious
    g_intsum(&total_mono_p[dir1]);
    g_intsum(&total_mono_m[dir1]);
  }

  total = 0;
  total_abs = 0;
  node0_printf("MONOPOLE ");
  for (dir1 = XUP; dir1 < NUMLINK - 1; dir1++) {
    total += total_mono_p[dir1] + total_mono_m[dir1];
    total_abs += total_mono_p[dir1] - total_mono_m[dir1];
    node0_printf("%d %d  ", total_mono_p[dir1], total_mono_m[dir1]);
  }
  node0_printf("  %d %d\n", total, total_abs);

  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    free(charge[dir1]);
    free(phase[dir1]);
    for (dir2 = 0; dir2 < NUMLINK; dir2++)
      free(mono[dir1][dir2]);
  }
}
#endif
// -----------------------------------------------------------------
