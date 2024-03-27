// -----------------------------------------------------------------
// Measure density of U(1) monopole world lines in doubled 3-cells
#include "susy_includes.h"

void monopole() {
  register int i, dir, dir2;
  register site *s;
  int a, b, c, d, ip2;
  int total_mono_p, total_mono_m, total = 0, total_abs = 0;
  int *mono[NDIMS][NDIMS], *charge;
  Real *phase[NDIMS], p2, p3, total_phase, permm;
  double threePI = 3.0 * PI, fivePI = 5.0 * PI, sevenPI = 7.0 * PI;
  complex det_link;
  msg_tag *mtag0, *mtag1;

  charge = malloc(sizeof(int) * sites_on_node);

  FORALLUPDIR(dir) {
    phase[dir] = malloc(sizeof(Real) * sites_on_node);
  }
  FORALLUPDIR(dir) {
    FORALLUPDIR(dir2) {
      mono[dir][dir2] = malloc(sizeof(int) * sites_on_node);
      FORALLSITES(i, s)
        mono[dir][dir2][i] = 0;     // Make sure we're initialized
    }
  }

  // First extract the U(1) part of the link
  FORALLUPDIR(dir) {
    FORALLSITES(i, s) {
      det_link = find_det(&(s->link[dir]));
      phase[dir][i] = atan2(det_link.imag, det_link.real);
//      printf("XXX (%d, %d, %d, %d)[%d] %.6g %.6g %.6g\n",
//             s->x, s->y, s->z, s->t, dir,
//             det_link.real, det_link.imag, phase[dir][i]);
    }
  }

  // Next find the number of strings
  // penetrating the face of every plaquette
  for (dir = YUP; dir <= TUP; dir++) {    // Omit fifth link
    for (dir2 = XUP; dir2 < dir; dir2++) {
      mtag0 = start_gather_field(phase[dir2], sizeof(Real),
                                 goffset[dir], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_field(phase[dir], sizeof(Real),
                                 goffset[dir2], EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        p2 = *((Real *)(gen_pt[0][i]));
        p3 = *((Real *)(gen_pt[1][i]));
        total_phase = phase[dir][i] + p2 - p3 - phase[dir2][i];
//        printf("YYY (%d, %d, %d, %d)[%d,%d] %.6g %.6g %.6g %.6g %.6g\n",
//               s->x, s->y, s->z, s->t, dir, dir2,
//               phase[dir][i], phase[dir2][i], p2, p3, total_phase);

        if (fabs(total_phase) < PI)
          mono[dir][dir2][i] = 0;
        else if (total_phase >= PI && total_phase < threePI)
          mono[dir][dir2][i] = 1;
        else if (total_phase >= threePI && total_phase < fivePI)
          mono[dir][dir2][i] = 2;
        else if (total_phase >= fivePI && total_phase < sevenPI)
          mono[dir][dir2][i] = 3;
        else if (total_phase <= -PI && total_phase > -threePI)
          mono[dir][dir2][i] = -1;
        else if (total_phase <= -threePI && total_phase > -fivePI)
          mono[dir][dir2][i] = -2;
        else if (total_phase <= -fivePI && total_phase > -sevenPI)
          mono[dir][dir2][i] = -3;
        else {
          printf("monopole: total_phase %.4g out of bounds on node%d\n",
                 total_phase, this_node);
          terminate(1);
        }
        mono[dir2][dir][i] = -mono[dir][dir2][i];

//        if (s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0) {
//          printf("ZZZ[%d, %d] ", dir, dir2);
//          printf("%.6g - %.6g + %.6g - %.6g = %.6g --> %d\n",
//                 phase[dir][i], phase[dir2][i], p2, p3,
//                 total_phase, mono[dir][dir2][i]);
//        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }

  // We have the number of strings penetrating every plaquette
  // Now tie these together into cubes,
  FORALLSITES(i, s)
    charge[i] = 0;

  FORALLUPDIR(a) {
    FORALLUPDIR(b) {
      if (b == a)
        continue;
      FORALLUPDIR(c) {
        if (c == a || c == b)
          continue;

        permm = perm[a][b][c];
//          printf("CHECK perm[%d][%d][%d] = %.1g \n",
//                 a, b, c, permm);

        mtag0 = start_gather_field(mono[b][c], sizeof(int),
                                   goffset[a], EVENANDODD, gen_pt[0]);

        wait_gather(mtag0);
        FORALLSITES(i, s) {
          ip2 = *((int *)(gen_pt[0][i]));
          if (permm > 0.0)
            charge[i] += (mono[b][c][i] - ip2);
          else if (permm < 0.0)
            charge[i] -= (mono[b][c][i] - ip2);
//            if (s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0) {
//              printf("charge[%d, %d, %d, %d] + (%d - %d) * (%d) --> %d\n",
//                     a, b, c, d, mono[c][d][i], ip2,
//                     (int)permm, charge[a][i]);
//            }
        }
        cleanup_gather(mtag0);
      }
    }
  }
  FORALLSITES(i, s)       // Normalize for epsilon loops double-counting
    charge[i] *= 0.5;  // Should end up an integer
//  FORALLSITES(i, s) {
//    FORALLUPDIR(dir) {
//      if (charge[dir][i] != 0) {
//        printf("QQQ (%d, %d, %d, %d)[%d] %d\n",
//               s->x, s->y, s->z, s->t, dir, charge[dir][i]);
//      }
//    }
//  }

  // Finally accumulate and print global quantities
  total_mono_p = 0;
  total_mono_m = 0;
  FORALLSITES(i, s) {
    if (charge[i] > 0)
      total_mono_p += charge[i];
    if (charge[i] < 0)
      total_mono_m += charge[i];
  }
  g_intsum(&total_mono_p);
  g_intsum(&total_mono_m);

  node0_printf("MONOPOLE ");
  
  if (total_mono_p + total_mono_m != 0)
    node0_printf("\nWARNING: total_mono mismatch in dir \n");
  total += total_mono_p + total_mono_m;
  total_abs += total_mono_p - total_mono_m;
  node0_printf("%d %d  ", total_mono_p, total_mono_m);

  node0_printf("  %d %d\n", total, total_abs);

  free(charge);
  FORALLUPDIR(dir) {
    free(phase[dir]);
    FORALLUPDIR(dir2)
      free(mono[dir][dir2]);
  }
}
// -----------------------------------------------------------------