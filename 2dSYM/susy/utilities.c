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
  int mu, nu, j, k;
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
        sub_su3_matrix_f(&(s->tempmat1), &(s->tempmat2), &(s->Fmunu[mu][nu]));
        for (j = 0; j < NCOL; j++) {
          for (k = 0; k < NCOL; k++)
            CNEGATE(s->Fmunu[mu][nu].e[j][k], s->Fmunu[nu][mu].e[j][k]);
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
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
  Dplus(link_src, plaq_dest);   // Overwrites plaq_dest for mu != nu
  Dminus(plaq_src, link_dest);  // Overwrites link_dest
#endif

#ifdef SV
  DbplusStoL(site_src, link_dest2);   // Overwrites link_dest2
  for (mu = 0; mu < NUMLINK; mu++) {
    FORALLSITES(i, s)
      scalar_mult_add_su3_vector(&(link_dest[mu][i]), &(link_dest2[mu][i]),
                                 0.5, &(link_dest[mu][i]));
  }

  DbminusLtoS(link_src, site_dest);   // Overwrites site_dest
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
  su3_vector vtmp1, vtmp2, vtmp3, vtmp4;
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
        mult_su3_mat_vec(&(s->link[mu]), (su3_vector *)(gen_pt[0][i]),
                         &vtmp1);
        scalar_mult_su3_vector(&vtmp1, s->bc[mu], &vtmp1);

        mult_su3_vec_mat(&(src[nu][i]), (su3_matrix *)(gen_pt[1][i]),
                         &vtmp2);
        mult_su3_mat_vec(&(s->link[nu]), (su3_vector *)(gen_pt[2][i]),
                         &vtmp3);
        scalar_mult_su3_vector(&vtmp3, s->bc[nu], &vtmp3);

        mult_su3_vec_mat(&(src[mu][i]), (su3_matrix *)(gen_pt[3][i]),
                         &vtmp4);
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
  int mu, nu, j;
  su3_vector tvec, vtmp1, vtmp2, vtmp3;
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
          for (j = 0; j < DIMF; j++)
            CNEGATE(src[i].c[j], tvec.c[j]);
        }
        else
          tvec = src[i];
        mult_su3_vec_mat(&tvec, &(s->link[mu]), &(tsite[0][i]));
      }
      mtag1 = start_gather_field(tsite[0], sizeof(su3_vector),
                                 goffset[mu] + 1, EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        if (mu > nu) {    // src is anti-symmetric under mu <--> nu
          for (j = 0; j < DIMF; j++)
            CNEGATE(src[i].c[j], tvec.c[j]);
        }
        else
          tvec = src[i];
        mult_su3_mat_vec((su3_matrix *)(gen_pt[0][i]), &tvec, &vtmp1);
        scalar_mult_su3_vector((su3_vector *)(gen_pt[1][i]),
                               s->bc[OPP_LDIR(mu)], &vtmp3);
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
  su3_vector vtmp1, vtmp2, vtmp3;
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
      mult_su3_vec_adj_mat(&(src[mu][i]), &(s->link[mu]), &vtmp1);
      scalar_mult_su3_vector((su3_vector *)(gen_pt[mu][i]),
                             s->bc[OPP_LDIR(mu)], &vtmp3);
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
  su3_vector vtmp1, vtmp2;
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
      mult_su3_vec_adj_mat((su3_vector *)(gen_pt[mu][i]), &(s->link[mu]),
                           &vtmp1);
      scalar_mult_su3_vector(&vtmp1, s->bc[mu], &vtmp1);
      mult_adj_su3_mat_vec(&(s->link[mu]), &(src[i]), &vtmp2);
      sub_su3_vector(&vtmp1, &vtmp2, &(dest[mu][i]));   // Overwrite
    }
    cleanup_gather(tag[mu]);
  }
}
// -----------------------------------------------------------------
