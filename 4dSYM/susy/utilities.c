// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force

// #define DET_DIST prints out all determinants for plotting distribution
// #define TR_DIST prints out all (Tr[U.Udag]/N-1)^2 for plotting distribution
// CAUTION: Do not run DET_DIST or TR_DIST with MPI!

//#define DET_DIST
//#define TR_DIST
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute at each site all NUMLINK * (NUMLINK - 1) plaquette determinants
// counting both orientations, and saving ZW* = plaqdet (plaqdet - 1)^*
// Use Tr_Uinv as temporary storage
void compute_plaqdet() {
  register int i;
  register site *s;
  int a, b;
  complex tc1, tc2;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  su3_matrix_f tmat;

#ifdef DET_DIST
  if (this_node != 0) {
    printf("compute_plaqdet: don't run DET_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  for (a = XUP; a < NUMLINK; a++) {
    for (b = XUP; b < NUMLINK; b++) {
      if (a == b)
        continue;

      // Gather determinants of U_a(x+b) and Udag_b(x+a)
      // as opposed to the full matrices
      // gen_pt[0] is det[Udag_b(x+a)], gen_pt[1] is det[U_a(x+b)]
      FORALLSITES(i, s) {
        Tr_Uinv[a][i] = find_det(&(s->linkf[a]));
        su3_adjoint_f(&(s->linkf[b]), &tmat);
        Tr_Uinv[b][i] = find_det(&tmat);
      }
      mtag0 = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                 goffset[a], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                 goffset[b], EVENANDODD, gen_pt[1]);

      // Compute det[U_b(x)] det[Udag_a(x)] while gathers run
      FORALLSITES(i, s) {
        su3_adjoint_f(&(s->linkf[a]), &tmat);
        tc1 = find_det(&tmat);
        tc2 = find_det(&(s->linkf[b]));
        CMUL(tc1, tc2, plaqdet[a][b][i]);       // Initialize
      }

      // Now put it all together
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        tc1 = *((complex *)(gen_pt[0][i]));     // det[Udag_b(x+a)]
        CMUL(plaqdet[a][b][i], tc1, tc2);
        tc1 = *((complex *)(gen_pt[1][i]));     // det[U_a(x+b)]
        CMUL(tc1, tc2, plaqdet[a][b][i]);

        // ZWstar = plaqdet (plaqdet - 1)^*
        CADD(plaqdet[a][b][i], minus1, tc1);
        CONJG(tc1, tc2);
        CMUL(plaqdet[a][b][i], tc2, ZWstar[a][b][i]);
#ifdef DET_DIST
        if (a < b) {
          printf("DET_DIST %d %d %d %d %d %d %.4g %.4g %.4g\n",
                 s->x, s->y, s->z, s->t, a, b,
                 plaqdet[a][b][i].real, plaqdet[a][b][i].imag, cabs_sq(&tc1));
        }
#endif
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Compute at each site
//   sum_mu [U_mu(x) * Udag_mu(x) - Udag_mu(x - mu) * U_mu(x - mu)]
// Add plaquette determinant contribution if G is non-zero
// Add scalar potential contribution if B is non-zero
// Use tempmat1 and tempmat2 as temporary storage
void compute_DmuUmu() {
  register int i;
  register site *s;
  int mu, nu, j;
  Real tr;
  complex tc;
  msg_tag *mtag0 = NULL;
  su3_matrix_f tmat, *mat;

#ifdef TR_DIST
  if (this_node != 0) {
    printf("compute_DmuUmu: don't run TR_DIST in parallel\n");
    fflush(stdout);
    terminate(1);
  }
#endif

  for (mu = XUP; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      mult_su3_na_f(&(s->linkf[mu]), &(s->linkf[mu]), &(tempmat1[i]));
      mult_su3_an_f(&(s->linkf[mu]), &(s->linkf[mu]), &(tempmat2[i]));
    }

    // Gather tempmat2 from below
    mtag0 = start_gather_field(tempmat2, sizeof(su3_matrix_f),
                               goffset[mu] + 1, EVENANDODD, gen_pt[0]);
    wait_gather(mtag0);
    if (mu == 0) {
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[0][i]);
        sub_su3_matrix_f(&(tempmat1[i]), mat, &(DmuUmu[i]));
      }
    }
    else {
      FORALLSITES(i, s) {
        mat = (su3_matrix_f *)(gen_pt[0][i]);
        sub_su3_matrix_f(&(tempmat1[i]), mat, &tmat);
        add_su3_matrix_f(&(DmuUmu[i]), &tmat, &(DmuUmu[i]));
      }
    }
    cleanup_gather(mtag0);
  }

  // Add plaquette determinant contribution if G is non-zero
  if (G > IMAG_TOL) {
    compute_plaqdet();
    FORALLSITES(i, s) {
      for (mu = XUP; mu < NUMLINK; mu++) {
        for (nu = XUP; nu < NUMLINK; nu++) {
          if (mu == nu)
            continue;

          CADD(plaqdet[mu][nu][i], minus1, tc);
#ifdef LINEAR_DET
          CMULREAL(tc, G, tc);
          for (j = 0; j < NCOL; j++)
            CADD(DmuUmu[i].e[j][j], tc, DmuUmu[i].e[j][j]);
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
  if (B > IMAG_TOL) {
    FORALLSITES(i, s) {
      for (mu = XUP; mu < NUMLINK; mu++) {
        tr = 1.0 / (Real)NCOL;
        tr *= realtrace_su3_f(&(s->linkf[mu]), &(s->linkf[mu]));
        tr -= 1.0;
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
// Compute at each site
//   U_mu(x) * U_mu(x + mu) - Udag_nu(x) * U_mu(x + nu)
// Use tempmat1 and tempmat2 as temporary storage
void compute_Fmunu() {
  register int i;
  register site *s;
  int mu, nu, index;
  su3_matrix_f *mat0, *mat1;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;

  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      index = plaq_index[mu][nu];
      mtag0 = start_gather_site(F_OFFSET(linkf[nu]), sizeof(su3_matrix_f),
                                goffset[mu], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[mu]), sizeof(su3_matrix_f),
                                goffset[nu], EVENANDODD, gen_pt[1]);
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mat0 = (su3_matrix_f *)(gen_pt[0][i]);
        mat1 = (su3_matrix_f *)(gen_pt[1][i]);
        mult_su3_nn_f(&(s->linkf[mu]), mat0, &(tempmat1[i]));
        mult_su3_nn_f(&(s->linkf[nu]), mat1, &(tempmat2[i]));
        sub_su3_matrix_f(&(tempmat1[i]), &(tempmat2[i]), &(Fmunu[index][i]));
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute traces of bilinears of scalar field interpolating ops:
//   traceless part of the hermitian matrix returned by the polar projection
//   traceless part of U_a Udag_a
#ifdef CORR
void compute_Ba() {
  register int i;
  register site *s;
  int a, b, j;
  complex tc;
  Real tr;
  su3_matrix_f tmat;

  // B_a is traceless part of the hermitian matrix from the polar projection
  // C_a is traceless part of U_a Udag_a
  FORALLSITES(i, s) {
    for (a = XUP; a < NUMLINK; a++) {
      polar(&(s->linkf[a]), &tmat, &(Ba[a][i]));
      tc = trace_su3_f(&(Ba[a][i]));
      tr = tc.real / (Real)NCOL;
      for (j = 0; j < NCOL; j++)
        Ba[a][i].e[j][j].real -= tr;

      mult_su3_na_f(&(s->linkf[a]), &(s->linkf[a]), &(Ca[a][i]));
      tc = trace_su3_f(&(Ca[a][i]));
      tr = tc.real / (Real)NCOL;
      for (j = 0; j < NCOL; j++)
        Ca[a][i].e[j][j].real -= tr;
    }

    // Bilinear traces are symmetric in a <--> b
    // traceBB[0][a][b] = tr[B_a(x) B_b(x)], etc.
    // But store all to simplify SUGRA computation
    // Make sure all are purely real
    for (a = XUP; a < NUMLINK; a++) {
      for (b = a; b < NUMLINK; b++) {
        mult_su3_nn_f(&(Ba[a][i]), &(Ba[b][i]), &tmat);
        tc = trace_su3_f(&tmat);
        traceBB[0][a][b][i] = tc.real;
        traceBB[0][a][b][i] = traceBB[0][b][a][i];
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr(BB[0][%d][%d]) = (%.4g, %.4g) at site %d\n",
                 a, b, tc.real, tc.imag, i);
        }

        mult_su3_nn_f(&(Ca[a][i]), &(Ca[b][i]), &tmat);
        tc = trace_su3_f(&tmat);
        traceBB[1][a][b][i] = tc.real;
        traceBB[1][a][b][i] = traceBB[1][b][a][i];
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr(BB[2][%d][%d]) = (%.4g, %.4g) at site %d\n",
                 a, b, tc.real, tc.imag, i);
        }

        // Cross term may not be so useful
        mult_su3_nn_f(&(Ba[a][i]), &(Ca[b][i]), &tmat);
        tc = trace_su3_f(&tmat);
        traceBB[2][a][b][i] = tc.real;
        traceBB[2][a][b][i] = traceBB[2][b][a][i];
        if (fabs(tc.imag) > IMAG_TOL) {
          printf("WARNING: Tr(BB[1][%d][%d]) = (%.4g, %.4g) at site %d\n",
                 a, b, tc.real, tc.imag, i);
        }
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate routines for each term in the fermion operator
// All called by fermion_op at the bottom of the file
#ifdef VP
void Dplus(su3_vector *src[NUMLINK], su3_vector *dest[NPLAQ]) {
  register int i;
  register site *s;
  int mu, nu, index;
  su3_vector vtmp1, vtmp2, vtmp3, vtmp4, *vec0, *vec2;
  su3_matrix *mat1, *mat3;
  msg_tag *mtag0 = NULL, *mtag1 = NULL, *mtag2 = NULL, *mtag3 = NULL;

  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      index = plaq_index[mu][nu];
      mtag0 = start_gather_field(src[nu], sizeof(su3_vector),
                                 goffset[mu], EVENANDODD, gen_pt[0]);

      mtag1 = start_gather_site(F_OFFSET(link[mu]), sizeof(su3_matrix),
                                goffset[nu], EVENANDODD, gen_pt[1]);

      mtag2 = start_gather_field(src[mu], sizeof(su3_vector),
                                 goffset[nu], EVENANDODD, gen_pt[2]);

      mtag3 = start_gather_site(F_OFFSET(link[nu]), sizeof(su3_matrix),
                                goffset[mu], EVENANDODD, gen_pt[3]);

      wait_gather(mtag0);
      wait_gather(mtag1);
      wait_gather(mtag2);
      wait_gather(mtag3);
      FORALLSITES(i, s) {
        vec0 = (su3_vector *)(gen_pt[0][i]);
        mat1 = (su3_matrix *)(gen_pt[1][i]);
        vec2 = (su3_vector *)(gen_pt[2][i]);
        mat3 = (su3_matrix *)(gen_pt[3][i]);
        mult_su3_mat_vec(&(s->link[mu]), vec0, &vtmp1);
        scalar_mult_su3_vector(&vtmp1, s->bc1[mu], &vtmp1);

        mult_su3_vec_mat(&(src[nu][i]), mat1, &vtmp2);
        mult_su3_mat_vec(&(s->link[nu]), vec2, &vtmp3);
        scalar_mult_su3_vector(&vtmp3, s->bc1[nu], &vtmp3);

        mult_su3_vec_mat(&(src[mu][i]), mat3, &vtmp4);
        sub_su3_vector(&vtmp1, &vtmp2, &vtmp1);
        sub_su3_vector(&vtmp3, &vtmp4, &vtmp3);
        sub_su3_vector(&vtmp1, &vtmp3, &(dest[index][i]));   // Initialize
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
      cleanup_gather(mtag2);
      cleanup_gather(mtag3);
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempvec[0] for temporary storage
#ifdef VP
void Dminus(su3_vector *src[NPLAQ], su3_vector *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu, nu, index;
  su3_vector vtmp1, vtmp2, *vec;
  su3_matrix *mat;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;

  for (nu = XUP; nu < NUMLINK; nu++) {
    FORALLSITES(i, s)
      clearvec(&(dest[nu][i]));         // Initialize

    for (mu = XUP; mu < NUMLINK; mu++) {
      if (mu == nu)
        continue;

      index = plaq_index[mu][nu];
      mtag0 = start_gather_site(F_OFFSET(link[mu]), sizeof(su3_matrix),
                                goffset[nu], EVENANDODD, gen_pt[0]);

      FORALLSITES(i, s) {
        if (mu > nu) {    // src is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(src[index][i]), -1.0, &vtmp1);
          mult_su3_vec_mat(&vtmp1, &(s->link[mu]), &(tempvec[0][i]));
        }
        else
          mult_su3_vec_mat(&(src[index][i]), &(s->link[mu]), &(tempvec[0][i]));
      }
      mtag1 = start_gather_field(tempvec[0], sizeof(su3_vector),
                                 goffset[mu] + 1, EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mat = (su3_matrix *)(gen_pt[0][i]);
        vec = (su3_vector *)(gen_pt[1][i]);
        if (mu > nu) {    // src is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(src[index][i]), -1.0, &vtmp2);
          mult_su3_mat_vec(mat, &vtmp2, &vtmp1);
        }
        else
          mult_su3_mat_vec(mat, &(src[index][i]), &vtmp1);
        scalar_mult_su3_vector(vec, s->bc1[OPP_LDIR(mu)], &vtmp2);
        sub_su3_vector(&vtmp1, &vtmp2, &vtmp1);
        add_su3_vector(&(dest[nu][i]), &vtmp1, &(dest[nu][i]));
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add to dest instead of overwriting; note factor of 1/2
#ifdef QCLOSED
void DbplusPtoP(su3_vector *src[NPLAQ], su3_vector *dest[NPLAQ]) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int a, b, c, d, e, j, gather, next, flip = 0, i_ab, i_de;
  Real tr;
  su3_vector *vec1, *vec2, vtmp1, vtmp2;
  su3_matrix *mat0, *mat3;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // Start first set of gathers
  // From setup_lambda.c, we see b > a and e > d
  c = DbplusPtoP_lookup[0][2];
  d = DbplusPtoP_lookup[0][3];
  e = DbplusPtoP_lookup[0][4];
  i_de = plaq_index[d][e];
  tag0[0] = start_gather_site(F_OFFSET(link[c]), sizeof(su3_matrix),
                              DbpP_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(src[i_de], sizeof(su3_vector),
                               DbpP_d2[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(src[i_de], sizeof(su3_vector),
                               DbpP_d1[0], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(link[c]), sizeof(su3_matrix),
                              goffset[c] + 1, EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {               // Start next set of gathers
      next = j + 1;
      c = DbplusPtoP_lookup[next][2];
      d = DbplusPtoP_lookup[next][3];
      e = DbplusPtoP_lookup[next][4];
      i_de = plaq_index[d][e];

      tag0[gather] = start_gather_site(F_OFFSET(link[c]),
                              sizeof(su3_matrix), DbpP_d1[next], EVENANDODD,
                              local_pt[gather][0]);
      tag1[gather] = start_gather_field(src[i_de],
                              sizeof(su3_vector), DbpP_d2[next], EVENANDODD,
                              local_pt[gather][1]);
      tag2[gather] = start_gather_field(src[i_de],
                              sizeof(su3_vector), DbpP_d1[next], EVENANDODD,
                              local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(link[c]),
                              sizeof(su3_matrix), goffset[c] + 1, EVENANDODD,
                              local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = DbplusPtoP_lookup[j][0];
    b = DbplusPtoP_lookup[j][1];
    c = DbplusPtoP_lookup[j][2];
    d = DbplusPtoP_lookup[j][3];
    e = DbplusPtoP_lookup[j][4];
    tr = 0.5 * perm[a][b][c][d][e];
    i_ab = plaq_index[a][b];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      mat0 = (su3_matrix *)(local_pt[flip][0][i]);
      vec1 = (su3_vector *)(local_pt[flip][1][i]);
      vec2 = (su3_vector *)(local_pt[flip][2][i]);
      mat3 = (su3_matrix *)(local_pt[flip][3][i]);
      mult_su3_vec_adj_mat(vec1, mat0, &vtmp1);
      scalar_mult_su3_vector(&vtmp1, s->bc3[a][b][c], &vtmp1);

      mult_adj_su3_mat_vec(mat3, vec2, &vtmp2);
      scalar_mult_su3_vector(&vtmp2, s->bc2[a][b], &vtmp2);

      sub_su3_vector(&vtmp1, &vtmp2, &vtmp1);
      scalar_mult_add_su3_vector(&(dest[i_ab][i]), &vtmp1, tr,
                                 &(dest[i_ab][i]));
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add to dest instead of overwriting; note factor of 1/2
#ifdef QCLOSED
void DbminusPtoP(su3_vector *src[NPLAQ], su3_vector *dest[NPLAQ]) {
  register int i;
  register site *s;
  char **local_pt[2][4];
  int a, b, c, d, e, j, gather, next, flip = 0, i_ab, i_de;
  Real tr;
  su3_vector *vec1, *vec2, vtmp1, vtmp2, vtmp3;
  su3_matrix *mat0, *mat3;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }

  // Start first set of gathers
  // From setup_lambda.c, we see b > a and e > d
  a = DbminusPtoP_lookup[0][0];
  b = DbminusPtoP_lookup[0][1];
  c = DbminusPtoP_lookup[0][2];
  i_ab = plaq_index[a][b];

  tag0[0] = start_gather_site(F_OFFSET(link[c]), sizeof(su3_matrix),
                              DbmP_d1[0], EVENANDODD, local_pt[0][0]);
  tag1[0] = start_gather_field(src[i_ab], sizeof(su3_vector),
                               DbmP_d2[0], EVENANDODD, local_pt[0][1]);
  tag2[0] = start_gather_field(src[i_ab], sizeof(su3_vector),
                               DbmP_d1[0], EVENANDODD, local_pt[0][2]);
  tag3[0] = start_gather_site(F_OFFSET(link[c]), sizeof(su3_matrix),
                              goffset[c] + 1, EVENANDODD, local_pt[0][3]);

  // Loop over lookup table
  for (j = 0; j < NTERMS; j++) {
    gather = (flip + 1) % 2;
    if (j < NTERMS - 1) {               // Start next set of gathers
      next = j + 1;
      a = DbminusPtoP_lookup[next][0];
      b = DbminusPtoP_lookup[next][1];
      c = DbminusPtoP_lookup[next][2];
      i_ab = plaq_index[a][b];
      tag0[gather] = start_gather_site(F_OFFSET(link[c]),
                              sizeof(su3_matrix), DbmP_d1[next], EVENANDODD,
                              local_pt[gather][0]);
      tag1[gather] = start_gather_field(src[i_ab],
                              sizeof(su3_vector), DbmP_d2[next], EVENANDODD,
                              local_pt[gather][1]);
      tag2[gather] = start_gather_field(src[i_ab],
                              sizeof(su3_vector), DbmP_d1[next], EVENANDODD,
                              local_pt[gather][2]);
      tag3[gather] = start_gather_site(F_OFFSET(link[c]),
                              sizeof(su3_matrix), goffset[c] + 1, EVENANDODD,
                              local_pt[gather][3]);
    }

    // Do this set of computations while next set of gathers runs
    a = DbminusPtoP_lookup[j][0];
    b = DbminusPtoP_lookup[j][1];
    c = DbminusPtoP_lookup[j][2];
    d = DbminusPtoP_lookup[j][3];
    e = DbminusPtoP_lookup[j][4];
    tr = 0.5 * perm[a][b][c][d][e];
    i_de = plaq_index[d][e];

    wait_gather(tag0[flip]);
    wait_gather(tag1[flip]);
    wait_gather(tag2[flip]);
    wait_gather(tag3[flip]);
    FORALLSITES(i, s) {
      mat0 = (su3_matrix *)(local_pt[flip][0][i]);
      vec1 = (su3_vector *)(local_pt[flip][1][i]);
      vec2 = (su3_vector *)(local_pt[flip][2][i]);
      mat3 = (su3_matrix *)(local_pt[flip][3][i]);
      mult_su3_vec_adj_mat(vec1, mat0, &vtmp1);
      scalar_mult_su3_vector(&vtmp1, s->bc2[OPP_LDIR(a)][OPP_LDIR(b)],
                             &vtmp1);

      mult_adj_su3_mat_vec(mat3, vec2, &vtmp2);
      scalar_mult_su3_vector(&vtmp2,
                             s->bc3[OPP_LDIR(a)][OPP_LDIR(b)][OPP_LDIR(c)],
                             &vtmp2);

      sub_su3_vector(&vtmp1, &vtmp2, &vtmp3);
      scalar_mult_add_su3_vector(&(dest[i_de][i]), &vtmp3, tr,
                                 &(dest[i_de][i]));
    }
    cleanup_gather(tag0[flip]);
    cleanup_gather(tag1[flip]);
    cleanup_gather(tag2[flip]);
    cleanup_gather(tag3[flip]);
    flip = (flip + 1) % 2;
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting site fermion to the link fermions
// bc1[mu](x) on psi_mu(x) eta(x + mu)
// Add to dest instead of overwriting; note factor of 1/2
#ifdef SV
void DbplusStoL(su3_vector *src, su3_vector *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  su3_vector vtmp1, vtmp2, *vec;
  msg_tag *tag[NUMLINK];

  tag[0] = start_gather_field(src, sizeof(su3_vector),
                              goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = XUP; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(su3_vector),
                                       goffset[mu + 1], EVENANDODD,
                                       gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      vec = (su3_vector *)(gen_pt[mu][i]);
      mult_su3_vec_adj_mat(vec, &(s->link[mu]), &vtmp1);
      scalar_mult_su3_vector(&vtmp1, s->bc1[mu], &vtmp1);
      mult_adj_su3_mat_vec(&(s->link[mu]), &(src[i]), &vtmp2);
      sub_su3_vector(&vtmp1, &vtmp2, &vtmp1);
      scalar_mult_add_su3_vector(&(dest[mu][i]), &vtmp1, 0.5, &(dest[mu][i]));
    }
    cleanup_gather(tag[mu]);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from site source to link destination
//   T[a](x) * sum_b {D[b][a](x) + D[a][b](x - b)}
// In the global case D is Tr[eta] * plaqdet,
// In the local case D is Tr[eta] plaqdet (plaqdet - 1)^*
// In both cases T is Tr[U^{-1} Lambda] and Tr[eta] = i sqrt(N) eta^D
// Use tr_dest and Tr_Uinv[0] for temporary storage
// bc1[b](x - b) = bc1[-b](x) on eta(x - b) psi_a(x)
// Add negative to dest instead of overwriting
// Negative sign is due to anti-commuting eta past psi
#ifdef SV
void detStoL(su3_vector *src, su3_vector *dest[NUMLINK]) {
  register int i;
  register site *s;
  int a, b, j, next;
  complex tc, tc2, tc3;
  complex Gc = cmplx(0.0, -1.0 * C2 * G * sqrt((Real)NCOL));
#ifdef LINEAR_DET
  CMULREAL(Gc, 0.5, Gc);                  // Since not squared
#endif
  su3_matrix_f tmat, tmat2;
  msg_tag *tag[NUMLINK];

  // compute_plaqdet saves plaqdet (plaqdet - 1)^* in ZWstar
  compute_plaqdet();

  // Save Tr[eta(x)] plaqdet[a][b](x) in modified plaqdet[a][b]
  // Save Tr[eta(x)] ZWstar[a][b](x) in modified plaqdet[a][b]
  for (a = XUP; a < NUMLINK; a++) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLSITES(i, s) {
#ifdef LINEAR_DET
        CMUL(src[i].c[DIMF - 1], plaqdet[a][b][i], tc);
#else
        CMUL(src[i].c[DIMF - 1], ZWstar[a][b][i], tc);
#endif
        plaqdet[a][b][i] = tc;

#ifdef LINEAR_DET
        CMUL(src[i].c[DIMF - 1], plaqdet[b][a][i], tc);
#else
        CMUL(src[i].c[DIMF - 1], ZWstar[b][a][i], tc);
#endif
        plaqdet[b][a][i] = tc;
      }
    }
  }

  // Now we gather modified plaqdet in both cases
  // Start first gather for (a, b) = (0, 1)
  tag[1] = start_gather_field(plaqdet[0][1], sizeof(complex),
                              goffset[1] + 1, EVENANDODD, gen_pt[1]);

  for (a = XUP; a < NUMLINK; a++) {
    // Initialize accumulator for sum over b
    FORALLSITES(i, s)
      tr_dest[i] = cmplx(0.0, 0.0);

    for (b = XUP; b < NUMLINK; b++) {
      if (a == b)
        continue;

      // Start next gather unless we're doing the last (a=4, b=3)
      next = b + 1;
      if (next < NUMLINK && a + b < 2 * NUMLINK - 3) {
        if (next == a)              // Next gather is actually (a, b + 2)
          next++;

        tag[next] = start_gather_field(plaqdet[a][next], sizeof(complex),
                                       goffset[next] + 1, EVENANDODD,
                                       gen_pt[next]);
      }
      else if (next == NUMLINK) {   // Start next gather (a + 1, 0)
        tag[0] = start_gather_field(plaqdet[a + 1][0], sizeof(complex),
                                    goffset[0] + 1, EVENANDODD, gen_pt[0]);
      }

      // Accumulate modified plaqdet[b][a](x) + plaqdet[a][b](x - b)
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        CMULREAL(tc, s->bc1[OPP_LDIR(b)], tc);
        CADD(plaqdet[b][a][i], tc, tc2);
        CSUM(tr_dest[i], tc2);
      }
      cleanup_gather(tag[b]);
    }

    // Compute Tr[U_a^{-1} Lambda^j] times sum
    FORALLSITES(i, s) {
      invert(&(s->linkf[a]), &tmat);
      CMUL(tr_dest[i], Gc, tc);
      for (j = 0; j < DIMF; j++) {
        mult_su3_nn_f(&tmat, &(Lambda[j]), &tmat2);
        tc2 = trace_su3_f(&tmat2);
        CMUL(tc, tc2, tc3);
        CSUM(dest[a][i].c[j], tc3);
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar potential coupling from site source to link destination
//   Udag[a](x) eta(x) (Tr[U_a(x) Udag_a(x)] / N - 1)
// Add negative to dest instead of overwriting
// Negative sign is due to anti-commuting eta past psi
#ifdef SV
void potStoL(su3_vector *src, su3_vector *dest[NUMLINK]) {
  register int i, a, j;
  register site *s;
  Real tr;
  complex tc, tc2, tc3;
  complex Bc = cmplx(0.0, -1.0 * C2 * B * B / sqrt((Real)NCOL));
  su3_matrix_f tmat;

  FORALLSITES(i, s) {
    for (a = XUP; a < NUMLINK; a++) {
      tr = 1.0 / (Real)NCOL;
      tr *= realtrace_su3_f(&(s->linkf[a]), &(s->linkf[a]));
      tr -= 1.0;
      CMULREAL(Bc, tr, tc);
      CMUL(src[i].c[DIMF - 1], tc, tc2);
      for (j = 0; j < DIMF; j++) {
        mult_su3_na_f(&(Lambda[j]), &(s->linkf[a]), &tmat);
        tc = trace_su3_f(&tmat);
        CMUL(tc, tc2, tc3);
        CSUM(dest[a][i].c[j], tc3);
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting the link fermions to the site fermion
// Given src psi_a, dest is Dbar_a psi_a (Eq. 63 in the arXiv:1108.1503)
// Use tempvec for temporary storage
// bc1[OPP_LDIR(mu)](x) on eta(x - mu) psi_mu(x - mu)
// Initialize dest; note factor of 1/2
#ifdef SV
void DbminusLtoS(su3_vector *src[NUMLINK], su3_vector *dest) {
  register int i;
  register site *s;
  int mu;
  su3_vector vtmp1, vtmp2, vtmp3, *vec;
  msg_tag *tag[NUMLINK];

  FORALLSITES(i, s) {         // Set up first gather
    clearvec(&(dest[i]));     // Initialize
    mult_adj_su3_mat_vec(&(s->link[0]), &(src[0][i]), &(tempvec[0][i]));
  }
  tag[0] = start_gather_field(tempvec[0], sizeof(su3_vector),
      goffset[0] + 1, EVENANDODD, gen_pt[0]);

  for (mu = 0; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1) {   // Start next gather
      FORALLSITES(i, s) {
        mult_adj_su3_mat_vec(&(s->link[mu + 1]), &(src[mu + 1][i]),
                             &(tempvec[mu + 1][i]));
      }
      tag[mu + 1] = start_gather_field(tempvec[mu + 1], sizeof(su3_vector),
                                       goffset[mu + 1] + 1, EVENANDODD,
                                       gen_pt[mu + 1]);
    }

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      vec = (su3_vector *)(gen_pt[mu][i]);
      mult_su3_vec_adj_mat(&(src[mu][i]), &(s->link[mu]), &vtmp1);
      scalar_mult_su3_vector(vec, s->bc1[OPP_LDIR(mu)], &vtmp3);
      sub_su3_vector(&vtmp1, &vtmp3, &vtmp2);
      scalar_mult_add_su3_vector(&(dest[i]), &vtmp2, 0.5, &(dest[i]));
    }
    cleanup_gather(tag[mu]);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from link source to site destination
//   sum_{a, b} D[a][b](x) * {T[b](x) +  T[a](x + b)}
// In the global case D is plaqdet
// In the local case D is plaqdet (plaqdet - 1)^*
// In both cases T is Tr[U^{-1} psi]
// bc1[b](x) on eta(x) psi_a(x + b)
// Use Tr_Uinv for temporary storage
// Add to dest instead of overwriting
// Has same sign as DbminusLtoS
#ifdef SV
void detLtoS(su3_vector *src[NUMLINK], su3_vector *dest) {
  register int i;
  register site *s;
  int a, b, j, next;
  complex tc, tc2;
  complex Gc = cmplx(0.0, C2 * G * sqrt((Real)NCOL));
#ifdef LINEAR_DET
  CMULREAL(Gc, 0.5, Gc);                  // Since not squared
#endif
  su3_matrix_f tmat, tmat2;
  msg_tag *tag[NUMLINK];

  compute_plaqdet();

  // Prepare Tr[U_a^{-1} psi_a] = sum_j Tr[U_a^{-1} Lambda^j] psi_a^j
  // and save in Tr_Uinv[a]
  for (a = XUP; a < NUMLINK; a++) {
    FORALLSITES(i, s) {
      invert(&(s->linkf[a]), &tmat);
      Tr_Uinv[a][i] = cmplx(0.0, 0.0);              // Initialize
      for (j = 0; j < DIMF; j++) {
        mult_su3_nn_f(&tmat, &(Lambda[j]), &tmat2);
        tc = trace_su3_f(&tmat2);
        CMUL(tc, src[a][i].c[j], tc2);
        CSUM(Tr_Uinv[a][i], tc2);                   // Accumulate
      }
    }
  }

  // Start first gather of Tr[U_a^{-1} psi_a] from x + b for (0, 1)
  tag[1] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                              goffset[1], EVENANDODD, gen_pt[1]);

  for (a = XUP; a < NUMLINK; a++) {
    for (b = XUP; b < NUMLINK; b++) {
      if (a == b)
        continue;

      // Start next gather unless we're doing the last (a=4, b=3)
      next = b + 1;
      if (next < NUMLINK && a + b < 2 * NUMLINK - 3) {
        if (next == a)              // Next gather is actually (a, b + 2)
          next++;
        tag[next] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                       goffset[next], EVENANDODD,
                                       gen_pt[next]);
      }
      else if (next == NUMLINK) {   // Start next gather (a + 1, 0)
        tag[0] = start_gather_field(Tr_Uinv[a + 1], sizeof(complex),
                                    goffset[0], EVENANDODD, gen_pt[0]);
      }

      // Accumulate D[a][b](x) {T[b](x) + T[a](x + b)}
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        CMULREAL(tc, s->bc1[b], tc);
        CADD(Tr_Uinv[b][i], tc, tc2);
#ifdef LINEAR_DET
        CMUL(plaqdet[a][b][i], tc2, tc);
#else
        CMUL(ZWstar[a][b][i], tc2, tc);
#endif
        CMUL(tc, Gc, tc2);
        CSUM(dest[i].c[DIMF - 1], tc2);
      }
      cleanup_gather(tag[b]);
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Scalar potential coupling from link source to site destination
//   sum_a (Tr[U_a(x) Udag_a(x)] / N - 1)^2 psi(x) Udag[a](x)
// Add to dest instead of overwriting
// Has same sign as DbminusLtoS
#ifdef SV
void potLtoS(su3_vector *src[NUMLINK], su3_vector *dest) {
  register int i, a, j;
  register site *s;
  Real tr;
  complex tc, tc2, tc3;
  complex Bc = cmplx(0.0, C2 * B * B / sqrt((Real)NCOL));
  su3_matrix_f tmat;

  FORALLSITES(i, s) {
    for (a = XUP; a < NUMLINK; a++) {
      tr = 1.0 / (Real)NCOL;
      tr *= realtrace_su3_f(&(s->linkf[a]), &(s->linkf[a]));
      tr -= 1.0;
      CMULREAL(Bc, tr, tc2);
      for (j = 0; j < DIMF; j++) {
        mult_su3_na_f(&(Lambda[j]), &(s->linkf[a]), &tmat);
        tc = trace_su3_f(&tmat);
        CMUL(tc, src[a][i].c[j], tc3);
        CMUL(tc2, tc3, tc);
        CSUM(dest[i].c[DIMF - 1], tc);
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Twist_Fermion matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign) {
  register int i, j;
  register site *s;
  int mu;
  Twist_Fermion tf;

  // Copy src TwistFermion into fieldwise site, link and plaq fermions
  // All of the latter are overwritten -- don't need to clear explicitly
  if (sign == 1) {
    FORALLSITES(i, s) {
      for (mu = 0; mu < NPLAQ; mu++)
        plaq_src[mu][i] = src[i].Fplaq[mu];
      for (mu = XUP; mu < NUMLINK; mu++)
        link_src[mu][i] = src[i].Flink[mu];

      site_src[i] = src[i].Fsite;
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      conjTF(&(src[i]), &tf);
      for (mu = 0; mu < NPLAQ; mu++)
        plaq_src[mu][i] = tf.Fplaq[mu];
      for (mu = XUP; mu < NUMLINK; mu++)
        link_src[mu][i] = tf.Flink[mu];

      site_src[i] = tf.Fsite;
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }

  // Assemble separate routines for each term in the fermion operator
#ifdef VP
  Dplus(link_src, plaq_dest);             // Initializes plaq_dest
  Dminus(plaq_src, link_dest);            // Initializes link_dest
#endif

#ifdef SV
  DbplusStoL(site_src, link_dest);        // Adds to link_dest2

  // Site-to-link plaquette determinant contribution if G is non-zero
  if (G > IMAG_TOL)
    detStoL(site_src, link_dest);         // Adds to link_dest

  // Site-to-link scalar potential contribution if B is non-zero
  if (B > IMAG_TOL)
    potStoL(site_src, link_dest);         // Adds to link_dest

  DbminusLtoS(link_src, site_dest);       // Initializes site_dest

  // Link-to-site plaquette determinant contribution if G is non-zero
  if (G > IMAG_TOL)
    detLtoS(link_src, site_dest);         // Adds to site_dest

  // Link-to-site scalar potential contribution if B is non-zero
  if (B > IMAG_TOL)
    potLtoS(link_src, site_dest);         // Adds to site_dest
#endif

#ifdef QCLOSED
  DbminusPtoP(plaq_src, plaq_dest);       // Adds to plaq_dest
  DbplusPtoP(plaq_src, plaq_dest);        // Adds to plaq_dest
#endif

  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      for (mu = 0; mu < NPLAQ; mu++)
        dest[i].Fplaq[mu] = plaq_dest[mu][i];
      for (mu = XUP; mu < NUMLINK; mu++)
        dest[i].Flink[mu] = link_dest[mu][i];

      dest[i].Fsite = site_dest[i];
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      for (j = 0; j < DIMF; j++) {
        for (mu = 0; mu < NPLAQ; mu++)
          CNEGATE(plaq_dest[mu][i].c[j], tf.Fplaq[mu].c[j]);
        for (mu = XUP; mu < NUMLINK; mu++)
          CNEGATE(link_dest[mu][i].c[j], tf.Flink[mu].c[j]);
        CNEGATE(site_dest[i].c[j], tf.Fsite.c[j]);
      }
      conjTF(&tf, &(dest[i]));
    }
  }
}
// -----------------------------------------------------------------
