// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force

// #define DET_DIST prints out all determinants for plotting distribution
// CAUTION: Do not run DET_DIST with MPI!

//#define DET_DIST
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute at each site all NUMLINK * (NUMLINK - 1) plaquette determinants
// counting both orientations, and saving ZW* = plaqdet (plaqdet - 1)^*
// Use Tr_Uinv as temporary storage
void compute_plaqdet() {
  register int i;
  register site *s;
  char **local_pt[2][2];
  int a, b, flip = 0;
  complex tc;
  msg_tag *tag0[2], *tag1[2];

#ifdef DET_DIST
  if (this_node != 0) {
    printf("compute_plaqdet: don't run DET_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  for (a = 0; a < 2; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[2 + a];
  }

  // Gather determinants rather than the full matrices
  // Recall det[Udag] = (det[U])^*
  FORALLSITES(i, s) {
    FORALLDIR(a)
      Tr_Uinv[a][i] = find_det(&(s->link[a]));
  }

  // Start first set of gathers (a = 0 and b = 1)
  // local_pt[0][0] is det[U_1(x+0)], local_pt[0][1] is det[U_0(x+1)]
  tag0[0] = start_gather_field(Tr_Uinv[1], sizeof(complex),
                               goffset[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                               goffset[1], EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(a) {
    FORALLDIR(b) {
      if (a == b)
        continue;

      if (a == 0) {         // Start other set of gathers (a = 1 and b = 0)
        tag0[1] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                                     goffset[1], EVENANDODD, local_pt[1][0]);
        tag1[1] = start_gather_field(Tr_Uinv[1], sizeof(complex),
                                     goffset[0], EVENANDODD, local_pt[1][1]);
      }

      // Initialize plaqdet[a][b] with det[U_b(x)] det[Udag_a(x)]
      FORALLSITES(i, s)
        CMULJ_(Tr_Uinv[a][i], Tr_Uinv[b][i], plaqdet[a][b][i]);

      // Now put it all together
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        // local_pt[flip][0] is det[U_b(x+a)]
        // Conjugate it to get det[Udag_b(x+a)]
        CMUL_J(plaqdet[a][b][i], *((complex *)(local_pt[flip][0][i])), tc);
        // local_pt[flip][1] is det[U_a(x+b)]
        CMUL(*((complex *)(local_pt[flip][1][i])), tc, plaqdet[a][b][i]);

        // ZWstar = plaqdet (plaqdet - 1)^*
        CADD(plaqdet[a][b][i], minus1, tc);
        CMUL_J(plaqdet[a][b][i], tc, ZWstar[a][b][i]);
#ifdef DET_DIST
        if (a < b) {
          printf("DET_DIST %d %d %d %d %.4g %.4g %.4g\n",
                 s->x, s->t, a, b,
                 plaqdet[a][b][i].real, plaqdet[a][b][i].imag, cabs_sq(&tc1));
        }
#endif
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip++;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting site fermion to the link fermions
// bc[mu](x) on psi_mu(x) eta(x + mu)
// Add to dest instead of overwriting; note factor of 1/2
#ifdef SV
void DbplusStoL(matrix *src, matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;

  tag[0] = start_gather_field(src, sizeof(matrix), goffset[0],
                              EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu == 0)     // Start other gather
      tag[1] = start_gather_field(src, sizeof(matrix), goffset[1],
                                  EVENANDODD, gen_pt[1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_na((matrix *)(gen_pt[mu][i]), &(s->link[mu]), &tmat);
      scalar_mult_matrix(&tmat, s->bc[mu], &tmat);
      mult_an_dif(&(s->link[mu]), &(src[i]), &tmat);
      scalar_mult_sum_matrix(&tmat, 0.5, &(dest[mu][i]));
    }
    cleanup_gather(tag[mu]);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting the link fermions to the site fermion
// Given src psi_a, dest is Dbar_a psi_a (Eq. 63 in the arXiv:1108.1503)
// Use tempmat and tempmat2 for temporary storage
// bc[OPP_LDIR(mu)](x) on eta(x - mu) psi_mu(x - mu)
// Initialize dest; note factor of 1/2
#ifdef SV
void DbminusLtoS(matrix *src[NUMLINK], matrix *dest) {
  register int i, mu, opp_mu;
  register site *s;
  msg_tag *tag[NUMLINK];

  FORALLSITES(i, s) {         // Set up first gather
    clear_mat(&(dest[i]));     // Initialize
    mult_an(&(s->link[0]), &(src[0][i]), &(tempmat[i]));
  }
  tag[0] = start_gather_field(tempmat, sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);

  FORALLDIR(mu) {
    if (mu == 0) {                  // Start other gather
      FORALLSITES(i, s)
        mult_an(&(s->link[1]), &(src[1][i]), &(tempmat2[i]));
      tag[1] = start_gather_field(tempmat2, sizeof(matrix),
                                  goffset[1] + 1, EVENANDODD, gen_pt[1]);
    }

    opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]), s->bc[opp_mu],
                             &(dest[i]));
      mult_na_sum(&(src[mu][i]), &(s->link[mu]), &(dest[i]));
    }
    cleanup_gather(tag[mu]);
  }

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix(&(dest[i]), 0.5, &(dest[i]));
}
#endif
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
void fermion_op(matrix *src[NFERMION], matrix *dest[NFERMION], int sign) {
  register int i, mu;
  register site *s;

  // Copy src matrix into fieldwise site, link and plaq fermions,
  // overwriting all of the latter
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(src[i].Fsite), &(site_src[i]));
      FORALLDIR(mu)
        mat_copy(&(src[i].Flink[mu]), &(link_src[mu][i]));
      mat_copy(&(src[i].Fplaq), &(plaq_src[i]));
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      adjoint(&(src[i].Fsite), &(site_src[i]));
      FORALLDIR(mu)
        adjoint(&(src[i].Flink[mu]), &(link_src[mu][i]));
      adjoint(&(src[i].Fplaq), &(plaq_src[i]));
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }
  FORALLSITES(i, s)
    tr_eta[i] = trace(&(site_src[i]));

  // Assemble separate routines for each term in the fermion operator

#ifdef SV
  DbplusStoL(site_src, link_dest);        // Adds to link_dest

  // Site-to-link plaquette determinant contribution if G is non-zero
  // Only depends on Tr[eta(x)]

  DbminusLtoS(link_src, site_dest);       // Overwrites site_dest

#endif

  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(site_dest[i]), &(dest[i].Fsite));
      FORALLDIR(mu)
        mat_copy(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
      mat_copy(&(plaq_dest[i]), &(dest[i].Fplaq));
    }
  }
  else if (sign == -1) {    // Both negate and conjugate
    FORALLSITES(i, s) {
      neg_adjoint(&(site_dest[i]), &(dest[i].Fsite));
      FORALLDIR(mu)
        neg_adjoint(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
      neg_adjoint(&(plaq_dest[i]), &(dest[i].Fplaq));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Squared matrix--vector operation
//   dest = (D^2).src
// Use temp_ferm for temporary storage
void DSq(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  register int i;
  register site *s;

  fermion_op(src, temp_ferm, PLUS);
  fermion_op(temp_ferm, dest, MINUS);
}
// -----------------------------------------------------------------
