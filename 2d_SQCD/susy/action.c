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
// Use tempmat and tempmat2 as temporary storage
void compute_DmuUmu() {
  register int i;
  register site *s;
  int mu, nu, j;
  complex tc;
  msg_tag *mtag0 = NULL;

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
          CMULREAL(tc, G, tc);
          for (j = 0; j < NCOL; j++)
            CSUM(DmuUmu[i].e[j][j], tc);
        }
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
// For the gauge action, compute at each site
// phi(x)*phibar(x)-r*I
// Use tempmat for temporary storage
void compute_PhiSq() {
  register int i;
  register site *s;

  FORALLSITES(i, s) {
    fun_mult_na(&(s->funlink), &(s->funlink), &(PhiSq[i]));
    scalar_mult_dif_matrix(&(Identity), R, &(PhiSq[i]));
  }
}

// -----------------------------------------------------------------


// -----------------------------------------------------------------
// For the gauge action, compute at each site
// Dmu phi = U_mu(x) phu(x+mu) - phi(x) U_mu(x+z)
// where U_mu(x+z) = FIdentity for both mu values

void compute_DmuPhi() {
  register int i;
  register site *s;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  //TODO:fix it so that F_OFFSET could be used for funlink
  //mtag0 = start_gather_site(F_OFFSET(funlink), sizeof(funmatrix),
  //                          goffset[0], EVENANDODD, gen_pt[0]);
  //mtag1 = start_gather_site(F_OFFSET(funlink), sizeof(funmatrix),
  //                          goffset[2], EVENANDODD, gen_pt[1]);
  
  mtag0 = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                            goffset[1], EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                            goffset[3], EVENANDODD, gen_pt[1]);

  wait_gather(mtag0);
  wait_gather(mtag1);

  FORALLSITES(i, s){
    fun_mat_mult((matrix *)(gen_pt[0][i]), &(s->funlink), &(tempmat[i]));
    fun_mat_mult((matrix *)(gen_pt[1][i]), &(s->funlink), &(tempmat2[i]));
  }

  mtag0 = start_gather_field(tempmat, sizeof(funmatrix),
                             goffset[0], EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_field(tempmat, sizeof(funmatrix),
                             goffset[2], EVENANDODD, gen_pt[1]);


  wait_gather(mtag0);
  wait_gather(mtag1);

  FORALLSITES(i, s){
    //fun_mat_mult(&(s->link[0]), (matrix *)(gen_pt[0][i]), &(tempmat[i]));
    //fun_sub_matrix(&(tempmat[i]), &(s->funlink), DmuPhi[0][i]);

    //fun_mat_mult(&(s->link[1]), (matrix *)(gen_pt[1][i]), &(tempmat[i]));
    //fun_sub_matrix(&(tempmat[i]), &(s->funlink), DmuPhi[1][i]);

    fun_sub_matrix((funmatrix *)(gen_pt[0][i]), &(s->funlink), DmuPhi[0][i]);
    fun_sub_matrix((funmatrix *)(gen_pt[1][i]), &(s->funlink), DmuPhi[1][i]);
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
  matrix tmat, tmat2, tmat3;

  FORALLSITES(i, s) {
    // d^2 term normalized by C2 / 2
    // DmuUmu includes the plaquette determinant contribution if G is non-zero
    mult_nn(&(DmuUmu[i]), &(DmuUmu[i]), &tmat);
    scalar_mult_matrix(&tmat, norm, &tmat);

    // F^2 term
    mult_an(&(Fmunu[i]), &(Fmunu[i]), &tmat2);
    scalar_mult_sum_matrix(&tmat2, 2.0, &tmat);

#ifdef SPHI
    // SQCD part of the gauge action
    // not certain about signs, numeric factors
    // phi^4 term
    mult_nn(&(PhiSq[i]), &(PhiSq[i]), &tmat3);
    scalar_mult_sum_matrix(&tmat3, 0.5, &tmat);

    // DmuUmu phi^2 term
    mult_nn(&(DmuUmu[i]), &(PhiSq[i]), &tmat3);
    scalar_mult_dif_matrix(&tmat3, 1.0, &tmat);

    // DmuPhi^2 term
    // in 1505.00467 this term is given as
    // bar(DmuPhi)DmuPhi but this is a NCOLF NCOLF matrix
    // so here we calculate DmuPhi bar(DmuPhi)
    fun_mult_na(&(DmuPhi[0][i]), &(DmuPhi[0][i]), &tmat3);
    scalar_mult_sum_matrix(&tmat3, 1.0, &tmat);
    fun_mult_na(&(DmuPhi[1][i]), &(DmuPhi[1][i]), &tmat3);
    scalar_mult_sum_matrix(&tmat3, 1.0, &tmat);
#endif


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
  double sum = 0.0;
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
      sum += realtrace(&tmat, &tmat);
#else
      tr = one_ov_N * realtrace(&(s->link[a]), &(s->link[a])) - 1.0;
      sum += tr * tr;
#endif
    }
  }
  sum *= kappa * bmass * bmass;
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
// is 4 DIMF volume on average
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

  //TODO:talk about phi momenta
  FORALLSITES(i, s) {
    FORALLDIR(mu)
      sum += (double)realtrace(&(s->mom[mu]), &(s->mom[mu]));
    sum += (double)funa_realtrace(&(s->funmom), &(s->funmom));
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
