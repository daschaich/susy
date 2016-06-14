// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// For the gauge action and force, compute at each site
//   sum_mu [U_mu(x) * Udag_mu(x) - Udag_mu(x - mu) * U_mu(x - mu)]
// Add plaquette determinant contribution if G is non-zero
// Add scalar potential contribution if B is non-zero
// Use tempmat and tempmat2 as temporary storage
void compute_DmuUmu() {
  register int i;
  register site *s;
  int mu, nu, j;
  Real tr;
  complex tc;
  msg_tag *mtag0 = NULL;

#ifdef TR_DIST
  if (this_node != 0) {
    printf("compute_DmuUmu: don't run TR_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  FORALLDIR(mu) {
    FORALLSITES(i, s) {
      mult_na(&(s->link[mu]), &(s->link[mu]), &(tempmat[i]));
      mult_an(&(s->link[mu]), &(s->link[mu]), &(tempmat2[i]));
    }

    // Gather tempmat2 from below
    mtag0 = start_gather_field(tempmat2, sizeof(matrix),
                               goffset[mu] + 1, EVENANDODD, gen_pt[0]);
    wait_gather(mtag0);
    if (mu == 0) {
      FORALLSITES(i, s)        // Initialize
        sub_matrix(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &(DmuUmu[i]));
    }
    else {
      FORALLSITES(i, s) {
        sum_matrix(&(tempmat[i]), &(DmuUmu[i]));
        dif_matrix((matrix *)(gen_pt[0][i]), &(DmuUmu[i]));
      }
    }
    cleanup_gather(mtag0);
  }

  // Add plaquette determinant contribution if G is non-zero
  // Assume compute_plaqdet() has already been run
  if (doG) {
    FORALLSITES(i, s) {
      FORALLDIR(mu) {
        FORALLDIR(nu) {
          if (mu == nu)
            continue;

          CADD(plaqdet[mu][nu][i], minus1, tc);
#ifdef LINEAR_DET
          CMULREAL(tc, G, tc);
          for (j = 0; j < NCOL; j++)
            CSUM(DmuUmu[i].e[j][j], tc);
#else
          tr = G * cabs_sq(&tc);
          for (j = 0; j < NCOL; j++)
            DmuUmu[i].e[j][j].real += tr;
#endif
        }
      }
    }
  }

  // Add scalar potential contribution if B is non-zero
  if (doB) {
    FORALLSITES(i, s) {
      FORALLDIR(mu) {
        tr = one_ov_N * realtrace(&(s->link[mu]), &(s->link[mu])) - 1.0;
        for (j = 0; j < NCOL; j++)
          DmuUmu[i].e[j][j].real += B * B * tr * tr;
#ifdef TR_DIST
          printf("TR_DIST %d %d %d %d %d %.4g\n",
                 s->x, s->y, s->z, s->t, mu, tr);
#endif
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// For the gauge action and force, compute at each site
//   U_mu(x) * U_mu(x + mu) - Udag_nu(x) * U_mu(x + nu)
// Use tempmat and tempmat2 as temporary storage
void compute_Fmunu() {
  register int i;
  register site *s;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;

  mtag0 = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                            goffset[0], EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                            goffset[1], EVENANDODD, gen_pt[1]);
  wait_gather(mtag0);
  wait_gather(mtag1);
  FORALLSITES(i, s) {
    mult_nn(&(s->link[0]), (matrix *)(gen_pt[0][i]), &(tempmat[i]));
    mult_nn(&(s->link[1]), (matrix *)(gen_pt[1][i]), &(tempmat2[i]));
    sub_matrix(&(tempmat[i]), &(tempmat2[i]), &(Fmunu[i]));
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Standard gauge contribution to the action
// Include tunable coefficient C2 in the d^2 term
double gauge_action(int do_det) {
  register int i;
  register site *s;
  double g_action = 0.0, norm = 0.5 * C2;
  complex tc;
  matrix tmat, tmat2;

  FORALLSITES(i, s) {
    // d^2 term normalized by C2 / 2
    // DmuUmu includes the plaquette determinant contribution if G is non-zero
    // and the scalar potential contribution if B is non-zero
    mult_nn(&(DmuUmu[i]), &(DmuUmu[i]), &tmat);
    scalar_mult_matrix(&tmat, norm, &tmat);

    // F^2 term
    mult_an(&(Fmunu[i]), &(Fmunu[i]), &tmat2);
    scalar_mult_sum_matrix(&tmat2, 2.0, &tmat);

    if (do_det == 1) {
      det_project(&tmat, &tmat2);
      tc = trace(&tmat2);
    }
    else
      tc = trace(&tmat);

    g_action += tc.real;
#ifdef DEBUG_CHECK
    if (fabs(tc.imag) > IMAG_TOL)
      printf("node%d WARNING: Im[s_B[%d]] = %.4g\n", this_node, i, tc.imag);
#endif
  }
  g_action *= kappa;
  g_doublesum(&g_action);
  return g_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar potential contribution to the action
// Note factor of kappa
double bmass_action() {
  register int i, a;
  register site *s;
  double sum = 0.0, kappa_muSq = kappa * bmass * bmass;
#ifdef EIG_POT
  matrix tmat;
#else
  Real tr;
#endif

  FORALLSITES(i, s) {
    FORALLDIR(a) {
#ifdef EIG_POT
      mult_na(&(s->link[a]), &(s->link[a]), &tmat);
      scalar_add_diag(&tmat, -1.0);
      sum += kappa_muSq * realtrace(&tmat, &tmat);
#else
      tr = one_ov_N * realtrace(&(s->link[a]), &(s->link[a])) - 1.0;
      sum += kappa_muSq * tr * tr;
#endif
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant contribution to the action
// Assume compute_plaqdet() has already been run
double det_action() {
  register int i;
  register site *s;
  int a, b;
  double re, im, det_action = 0.0;

  FORALLSITES(i, s) {
    FORALLDIR(a) {
      for (b = a + 1; b < NUMLINK; b++) {
        re = plaqdet[a][b][i].real;
        im = plaqdet[a][b][i].imag;
        det_action += (re - 1.0) * (re - 1.0);
        det_action += im * im;
      }
    }
  }
  det_action *= kappa_u1;
  g_doublesum(&det_action);
  return det_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action
// Include the ampdeg term to allow sanity check that the fermion action
// is 16 DIMF volume on average
// Since the pseudofermion src is fixed throughout the trajectory,
// ampdeg actually has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double fermion_action(Twist_Fermion *src, Twist_Fermion **sol) {
  register int i, j;
  register site *s;
  double sum = 0.0;
  complex ctmp;
#ifdef DEBUG_CHECK
  double im = 0.0;
#endif

  FORALLSITES(i, s) {
    sum += ampdeg4 * (double)magsq_TF(&(src[i]));
    for (j = 0; j < Norder; j++) {
      ctmp = TF_dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (double)(amp4[j] * ctmp.real);
#ifdef DEBUG_CHECK  // Make sure imaginary part vanishes
      im += (double)(amp4[j] * ctmp.imag);
#endif
    }
  }
  g_doublesum(&sum);
#ifdef DEBUG_CHECK  // Make sure imaginary part vanishes
  g_doublesum(&im);
  node0_printf("S_f = (%.4g, %.4g)\n", sum, im);
#endif
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momenta contribution to the action
double mom_action() {
  register int i, mu;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s) {
    FORALLDIR(mu)
      sum += (double)realtrace(&(s->mom[mu]), &(s->mom[mu]));
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out zeros for pieces of the action that aren't included
double action(Twist_Fermion **src, Twist_Fermion ***sol) {
  double g_act, bmass_act = 0.0, p_act, f_act, det_act = 0.0;
  double total;

  g_act = gauge_action(NODET);
  if (bmass > IMAG_TOL)
    bmass_act = bmass_action();
  if (kappa_u1 > IMAG_TOL)
    det_act = det_action();

  node0_printf("action: ");
  node0_printf("gauge %.8g bmass %.8g ", g_act, bmass_act);
  node0_printf("det %.8g ", det_act);

  total = g_act + bmass_act + det_act;
#ifndef PUREGAUGE
  int n;
  for (n = 0; n < Nroot; n++) {
    f_act = fermion_action(src[n], sol[n]);
    node0_printf("fermion%d %.8g ", n, f_act);
    total += f_act;
  }
#endif
  p_act = mom_action();
  node0_printf("mom %.8g ", p_act);
  total += p_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
