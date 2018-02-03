// -----------------------------------------------------------------
// Measure density of U(1) monopole flux through plaquettes
#include "susy_includes.h"

void monopole() {
  register int i, dir;
  register site *s;
  int total_mono_p = 0, total_mono_m = 0, total = 0, total_abs = 0;
  int *mono = malloc(sizeof *mono * sites_on_node);
  Real *phase[NDIMS], p2, p3, total_phase;
  double threePI = 3.0 * PI, fivePI = 5.0 * PI, sevenPI = 7.0 * PI;
  complex det_link;
  msg_tag *mtag0, *mtag1;

  FORALLUPDIR(dir)
    phase[dir] = malloc(sizeof(Real) * sites_on_node);

  // First extract the U(1) part of the link
  FORALLUPDIR(dir) {
    FORALLSITES(i, s) {
      det_link = find_det(&(s->link[dir]));
      phase[dir][i] = atan2(det_link.imag, det_link.real);
//      printf("XXX (%d, %d)[%d] %.6g %.6g %.6g\n",
//             s->x, s->t, dir, det_link.real, det_link.imag, phase[dir][i]);
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

    if (fabs(total_phase) < PI)
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

  // Finally accumulate and print global quantities
  FORALLSITES(i, s) {
    if (mono[i] > 0)
      total_mono_p += mono[i];
    if (mono[i] < 0)
      total_mono_m += mono[i];
  }
  g_intsum(&total_mono_p);
  g_intsum(&total_mono_m);

  total += total_mono_p + total_mono_m;
  total_abs += total_mono_p - total_mono_m;
  if (total != 0)
    node0_printf("\nWARNING: total_mono mismatch\n");

  node0_printf("MONOPOLE %d %d  ", total_mono_p, total_mono_m);
  node0_printf("  %d %d\n", total, total_abs);

  FORALLUPDIR(dir)
    free(phase[dir]);
  free(mono);
}
// -----------------------------------------------------------------
