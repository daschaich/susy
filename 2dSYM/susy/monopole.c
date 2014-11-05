// -----------------------------------------------------------------
// Measure density of monopole world lines in non-diagonal cubes
#include "susy_includes.h"

void monopole() {
  register int i, dir, dir2;
  register site *s;
  int ip2, total_mono_p[NDIMS], total_mono_m[NDIMS], total, total_abs;
  int *charge[NDIMS];
  int *mono = malloc(sites_on_node * sizeof(*mono));
  Real *phase[NUMLINK], p2, p3, total_phase;
  double threePI = 3.0 * PI;
  double fivePI = 5.0 * PI;
  double sevenPI = 7.0 * PI;
  complex det_link;
  msg_tag *mtag0, *mtag1;

  for (dir = XUP; dir <= TUP; dir++) {
    charge[dir] = malloc(sites_on_node * sizeof(int));
    phase[dir] = malloc(sites_on_node * sizeof(Real));
  }

  // First extract the U(1) part of the link
  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, s) {
      det_link = find_det(&(s->linkf[dir]));
      phase[dir][i] = atan2(det_link.imag, det_link.real);
//      printf("XXX (%d, %d, %d, %d)[%d] %.6g %.6g %.6g\n",
//             s->x, s->y, s->z, s->t, dir,
//             det_link.real, det_link.imag, phase[dir][i]);
    }
  }

  // Next find the number of strings
  // penetrating the face of every plaquette
  // Only have one oriented plane
  mtag0 = start_gather_field(phase[XUP], sizeof(Real), goffset[TUP],
                             EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_field(phase[TUP], sizeof(Real), goffset[XUP],
                             EVENANDODD, gen_pt[1]);

  wait_gather(mtag0);
  wait_gather(mtag1);
  FORALLSITES(i, s) {
    p2 = *((Real *)(gen_pt[0][i]));
    p3 = *((Real *)(gen_pt[1][i]));
    total_phase = phase[TUP][i] + p2 - p3 - phase[XUP][i];
//    printf("YYY (%d, %d) %.6g %.6g %.6g %.6g %.6g\n",
//           s->x, s->t, phase[TUP][i], phase[XUP][i], p2, p3, total_phase);

    if (abs(total_phase) < PI)
      mono[i] = 0;
    else if (total_phase >= PI && total_phase < threePI)
      mono[i] = 1;
    else if (total_phase >= threePI && total_phase < fivePI)
      mono[i] = 2;
    else if (total_phase >= fivePI && total_phase < sevenPI)
      mono[i] = 3;
    else if (total_phase <= -PI && total_phase > -threePI)
      mono[i] = -1;
    else if (total_phase <= -threePI && total_phase > -fivePI)
      mono[i] = -2;
    else if (total_phase <= -fivePI && total_phase > -sevenPI)
      mono[i] = -3;
    else {
      printf("monopole: total_phase %.4g out of bounds on node%d\n",
             total_phase, this_node);
      terminate(1);
    }

//    if (s->x == 0 && s->t == 0) {
//      printf("ZZZ %.6g - %.6g + %.6g - %.6g = %.6g --> %d\n",
//             phase[TUP][i], phase[XUP][i], p2, p3, total_phase, mono[i]);
//    }
  }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);

  // We have the number of strings penetrating every plaquette
  // Now tie these together into squares
  for (dir = XUP; dir <= TUP; dir++) {
    dir2 = 1 - dir;
    FORALLSITES(i, s)
      charge[dir][i] = 0;   // Should end up an integer, potentially odd

    mtag0 = start_gather_field(mono, sizeof(int), goffset[dir2],
                               EVENANDODD, gen_pt[0]);
    wait_gather(mtag0);
    FORALLSITES(i, s) {
      ip2 = *((int *)(gen_pt[0][i]));
      if (dir2 > dir)
        charge[dir][i] += (mono[i] - ip2);
      else if (dir2 < dir)
        charge[dir][i] -= (mono[i] - ip2);
      else {
        printf("monopole: something weird happened\n");
        terminate(1);
      }
//      if (s->x == 0 && s->t == 0)
//        printf("%d %d %d %d\n", dir, dir2, mono[i], ip2);
    }
    cleanup_gather(mtag0);
  }
//  FORALLSITES(i, s) {
//    for (dir = XUP; dir < NUMLINK - 1; dir++) {
//      if (charge[dir][i] != 0) {
//        printf("QQQ (%d, %d, %d, %d)[%d] %d\n",
//               s->x, s->y, s->z, s->t, dir, charge[dir][i]);
//      }
//    }
//  }

  // Finally accumulate and print global quantities
  total = 0;
  total_abs = 0;
  for (dir = XUP; dir <= TUP; dir++) {
    total_mono_p[dir] = 0;
    total_mono_m[dir] = 0;
    FORALLSITES(i, s) {
      if (charge[dir][i] > 0)
        total_mono_p[dir] += charge[dir][i];
      if (charge[dir][i] < 0)
        total_mono_m[dir] += charge[dir][i];
    }
    g_intsum(&total_mono_p[dir]);
    g_intsum(&total_mono_m[dir]);
  }

  total = 0;
  total_abs = 0;
  node0_printf("MONOPOLE ");
  for (dir = XUP; dir <= TUP; dir++) {
    if (total_mono_p[dir] + total_mono_m[dir] != 0)
      node0_printf("\nWARNING: total_mono mismatch in dir %d\n", dir);
    total += total_mono_p[dir] + total_mono_m[dir];
    total_abs += total_mono_p[dir] - total_mono_m[dir];
    node0_printf("%d %d  ", total_mono_p[dir], total_mono_m[dir]);
  }
  node0_printf("  %d %d\n", total, total_abs);

  for (dir = XUP; dir <= TUP; dir++) {
    free(charge[dir]);
    free(phase[dir]);
  }
  free(mono);
}
// -----------------------------------------------------------------
