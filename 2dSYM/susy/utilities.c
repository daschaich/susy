// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute at each site
//   sum_mu [U_mu(x) * Udag_mu(x) - Udag_mu(x - mu) * U_mu(x - mu)]
// Use s->tempmat1, s->tempmat2 and s->staple as temporary storage
void compute_DmuUmu() {
  register int i;
  register site *s;
  int mu;
  msg_tag *mtag0 = NULL;

  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s) {
      mult_su3_na_f(&(s->linkf[mu]), &(s->linkf[mu]), &(s->tempmat1));
      mult_su3_an_f(&(s->linkf[mu]), &(s->linkf[mu]), &(s->tempmat2));
    }

    // Gather tempmat2 from below
    mtag0 = start_gather_site(F_OFFSET(tempmat2), sizeof(su3_matrix_f),
                              goffset[mu] + 1, EVENANDODD, gen_pt[0]);
    wait_gather(mtag0);
    if (mu == 0) {
      FORALLSITES(i, s)
        sub_su3_matrix_f(&(s->tempmat1), (su3_matrix_f *)(gen_pt[0][i]),
                         &(s->DmuUmu));
    }
    else {
      FORALLSITES(i, s) {
        sub_su3_matrix_f(&(s->tempmat1), (su3_matrix_f *)(gen_pt[0][i]),
                         &(s->staple));
        add_su3_matrix_f(&(s->DmuUmu), &(s->staple), &(s->DmuUmu));
      }
    }
    cleanup_gather(mtag0);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute at each site
//   U_mu(x) * U_mu(x + mu) - Udag_nu(x) * U_mu(x + nu)
// Use s->tempmat1 and s->tempmat2 as temporary storage
void compute_Fmunu() {
  register int i;
  register site *s;
  int mu, nu;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;

  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      mtag0 = start_gather_site(F_OFFSET(linkf[nu]), sizeof(su3_matrix_f),
                                goffset[mu], EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[mu]), sizeof(su3_matrix_f),
                                goffset[nu], EVENANDODD, gen_pt[1]);
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->linkf[mu]), (su3_matrix_f *)(gen_pt[0][i]),
                      &(s->tempmat1));
        mult_su3_nn_f(&(s->linkf[nu]), (su3_matrix_f *)(gen_pt[1][i]),
                      &(s->tempmat2));
        sub_su3_matrix_f(&(s->tempmat1), &(s->tempmat2), &(s->Fmunu));
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute at each site B_mu(x) = U_mu(x) * Udag_mu(x) - trace
// as well as traceBB[mu][nu](x) = tr[B_mu(x) B_nu(x)] (should be real)
#ifdef CORR
void compute_Bmu() {
  register int i;
  register site *s;
  int mu, nu, j;
  complex ctmp;
  su3_matrix_f tmat;

  // B_mu = U_mu Udag_mu - trace
  FORALLSITES(i, s) {
    for (mu = 0; mu < NUMLINK; mu++) {
      mult_su3_na_f(&(s->linkf[mu]), &(s->linkf[mu]), &(s->B[mu]));
      ctmp = trace_su3_f(&(s->B[mu]));
      CDIVREAL(ctmp, (Real)NCOL, ctmp);
      for (j = 0; j < NCOL; j++)
        CDIF(s->B[mu].e[j][j], ctmp);
    }

    // traceBB[mu][nu] = tr[B_mu(x) B_nu(x)] is symmetric in mu <--> nu
    // But store all to simplify SUGRA computation
    // Make sure Tr(B_a * B_b) is purely real
    for (mu = 0; mu < NUMLINK; mu++) {
      for (nu = mu; nu < NUMLINK; nu++) {
        mult_su3_nn_f(&(s->B[mu]), &(s->B[nu]), &tmat);
        ctmp = trace_su3_f(&tmat);
        if (fabs(ctmp.imag) > IMAG_TOL) {
          printf("node%d WARNING: Tr(BB[%d][%d]) = (%.4g, %.4g) at site %d\n",
                 this_node, mu, nu, ctmp.real, ctmp.imag, i);
        }
        s->traceBB[mu][nu] = ctmp.real;
        s->traceBB[nu][mu] = ctmp.real;
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
      plaq_src[i] = src[i].Fplaq;
      for (mu = 0; mu < NUMLINK; mu++)
        link_src[mu][i] = src[i].Flink[mu];

      site_src[i] = src[i].Fsite;
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      conjTF(&(src[i]), &tf);
      plaq_src[i] = tf.Fplaq;
      for (mu = 0; mu < NUMLINK; mu++)
        link_src[mu][i] = tf.Flink[mu];

      site_src[i] = tf.Fsite;
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op\n");
    exit(1);
  }

  // Now the fun begins
#ifdef VP
  Dplus(link_src, plaq_dest);             // Overwrites plaq_dest
  Dminus(plaq_src, link_dest);            // Overwrites link_dest
#endif

#ifdef SV
  DbplusStoL(site_src, link_dest2);       // Overwrites link_dest2
  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      scalar_mult_add_su3_vector(&(link_dest[mu][i]), &(link_dest2[mu][i]),
                                 0.5, &(link_dest[mu][i]));
  }

  DbminusLtoS(link_src, site_dest);       // Overwrites site_dest
  FORALLSITES(i, s)
    scalar_mult_su3_vector(&(site_dest[i]), 0.5, &(site_dest[i]));
#endif

  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      dest[i].Fplaq = plaq_dest[i];
      for (mu = 0; mu < NUMLINK; mu++)
        dest[i].Flink[mu] = link_dest[mu][i];

      dest[i].Fsite = site_dest[i];
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      for (j = 0; j < DIMF; j++) {
        CNEGATE(plaq_dest[i].c[j], tf.Fplaq.c[j]);
        for (mu = 0; mu < NUMLINK; mu++)
          CNEGATE(link_dest[mu][i].c[j], tf.Flink[mu].c[j]);
        CNEGATE(site_dest[i].c[j], tf.Fsite.c[j]);
      }
      conjTF(&tf, &(dest[i]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void Dplus(su3_vector *src[NUMLINK], su3_vector *dest) {
  register int i;
  register site *s;
  int mu, nu;
  su3_vector vtmp1, vtmp2, vtmp3, vtmp4, *vec0, *vec2;
  su3_matrix *mat1, *mat3;
  msg_tag *mtag0 = NULL, *mtag1 = NULL, *mtag2 = NULL, *mtag3 = NULL;

  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
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
        scalar_mult_su3_vector(&vtmp1, s->bc[mu], &vtmp1);

        mult_su3_vec_mat(&(src[nu][i]), mat1, &vtmp2);
        mult_su3_mat_vec(&(s->link[nu]), vec2, &vtmp3);
        scalar_mult_su3_vector(&vtmp3, s->bc[nu], &vtmp3);

        mult_su3_vec_mat(&(src[mu][i]), mat3, &vtmp4);
        sub_su3_vector(&vtmp1, &vtmp2, &vtmp1);
        sub_su3_vector(&vtmp3, &vtmp4, &vtmp3);
        sub_su3_vector(&vtmp1, &vtmp3, &(dest[i]));   // Overwrite
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
      cleanup_gather(mtag2);
      cleanup_gather(mtag3);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void Dminus(su3_vector *src, su3_vector *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu, nu;
  su3_vector vtmp1, vtmp2, vtmp3, *vec1, tvec;
  su3_matrix *mat0;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;

  for (nu = 0; nu < NUMLINK; nu++) {
    FORALLSITES(i, s)
      clearvec(&(dest[nu][i]));   // Overwrite

    for (mu = 0; mu < NUMLINK; mu++) {
      if (mu == nu)
        continue;

      mtag0 = start_gather_site(F_OFFSET(link[mu]), sizeof(su3_matrix),
                                goffset[nu], EVENANDODD, gen_pt[0]);

      FORALLSITES(i, s) {
        if (mu > nu) {    // src is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(src[i]), -1.0, &tvec);
          mult_su3_vec_mat(&tvec, &(s->link[mu]), &(tsite[0][i]));
        }
        else
          mult_su3_vec_mat(&(src[i]), &(s->link[mu]), &(tsite[0][i]));
      }
      mtag1 = start_gather_field(tsite[0], sizeof(su3_vector),
                                 goffset[mu] + 1, EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mat0 = (su3_matrix *)(gen_pt[0][i]);
        vec1 = (su3_vector *)(gen_pt[1][i]);
        if (mu > nu) {    // src is anti-symmetric under mu <--> nu
          scalar_mult_su3_vector(&(src[i]), -1.0, &tvec);
          mult_su3_mat_vec(mat0, &tvec, &vtmp1);
        }
        else
          mult_su3_mat_vec(mat0, &(src[i]), &vtmp1);
        scalar_mult_su3_vector(vec1, s->bc[OPP_LDIR(mu)], &vtmp3);
        sub_su3_vector(&vtmp1, &vtmp3, &vtmp2);
        add_su3_vector(&(dest[nu][i]), &vtmp2, &(dest[nu][i]));
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Term in action connecting the link fermions to a site fermion
// Given src psi_a, dest is Dbar_a psi_a (Eq. 63 in the arXiv:1108.1503)
void DbminusLtoS(su3_vector *src[NUMLINK], su3_vector *dest) {
  register int i;
  register site *s;
  int mu;
  su3_vector vtmp1, vtmp2, vtmp3, *vec;
  msg_tag *tag[NUMLINK];

  FORALLSITES(i, s) {         // Set up first gather
    clearvec(&(dest[i]));     // Overwrite
    mult_adj_su3_mat_vec(&(s->link[0]), &(src[0][i]), &(tsite[0][i]));
  }
  tag[0] = start_gather_field(tsite[0], sizeof(su3_vector),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);

  for (mu = 0; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1) {   // Start next gather
      FORALLSITES(i, s)
        mult_adj_su3_mat_vec(&(s->link[mu + 1]), &(src[mu + 1][i]),
                             &(tsite[mu + 1][i]));

      tag[mu + 1] = start_gather_field(tsite[mu + 1], sizeof(su3_vector),
                                       goffset[mu + 1] + 1, EVENANDODD,
                                       gen_pt[mu + 1]);
    }

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      vec = (su3_vector *)(gen_pt[mu][i]);
      mult_su3_vec_adj_mat(&(src[mu][i]), &(s->link[mu]), &vtmp1);
      scalar_mult_su3_vector(vec, s->bc[OPP_LDIR(mu)], &vtmp3);
      sub_su3_vector(&vtmp1, &vtmp3, &vtmp2);
      add_su3_vector(&(dest[i]), &vtmp2, &(dest[i]));
    }
    cleanup_gather(tag[mu]);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void DbplusStoL(su3_vector *src, su3_vector *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  su3_vector vtmp1, vtmp2, *vec;
  msg_tag *tag[NUMLINK];

  tag[0] = start_gather_field(src, sizeof(su3_vector),
                              goffset[0], EVENANDODD, gen_pt[0]);
  for (mu = 0; mu < NUMLINK; mu++) {
    if (mu < NUMLINK - 1)  // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(su3_vector),
                                       goffset[mu + 1], EVENANDODD,
                                       gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      vec = (su3_vector *)(gen_pt[mu][i]);
      mult_su3_vec_adj_mat(vec, &(s->link[mu]), &vtmp1);
      scalar_mult_su3_vector(&vtmp1, s->bc[mu], &vtmp1);
      mult_adj_su3_mat_vec(&(s->link[mu]), &(src[i]), &vtmp2);
      sub_su3_vector(&vtmp1, &vtmp2, &(dest[mu][i]));   // Overwrite
    }
    cleanup_gather(tag[mu]);
  }
}
// -----------------------------------------------------------------
