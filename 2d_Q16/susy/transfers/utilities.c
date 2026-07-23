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
  int a, b, gather, flip = 0, mu, nu;
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

      gather = (flip + 1) % 2;
      if (a < NUMLINK - 1 || b < NUMLINK - 2) { // Start next set of gathers
        if (b == NUMLINK - 1) {
          mu = a + 1;
          nu = 0;
        }
        else if (b == a - 1) {
          mu = a;
          nu = b + 2;
        }
        else {
          mu = a;
          nu = b + 1;
        }
        tag0[gather] = start_gather_field(Tr_Uinv[nu], sizeof(complex),
                                          goffset[mu], EVENANDODD,
                                          local_pt[gather][0]);
        tag1[gather] = start_gather_field(Tr_Uinv[mu], sizeof(complex),
                                          goffset[nu], EVENANDODD,
                                          local_pt[gather][1]);
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
      flip = gather;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Save U_a^{-1} and Udag_a^{-1} = (U_a^{-1})^dag at each site
void compute_Uinv() {
  register int i, mu;
  register site *s;

  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      invert(&(s->link[mu]), &(Uinv[mu][i]));
      adjoint(&(Uinv[mu][i]), &(Udag_inv[mu][i]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Separate routines for each term in the fermion operator
// All called by fermion_op at the bottom of the file
#ifdef VP
void Dplus(matrix *src[NUMLINK], matrix *dest[NPLAQ]) {

//edited 
      int p, q, indx;
 for (p = 0; p < NUMLINK; p++) {
    plaq_index[p][p] = -1;                                  // i,i=-1  ,  00=11=22=-1
    for (q = p + 1; q < NUMLINK; q++) {
      indx = p * (NUMLINK - 1) - p * (p + 1) / 2 + q - 1;//              01=0,02=1,12=2 
      plaq_index[p][q] = indx;
      plaq_index[q][p] = indx;                             // i,j=j,i, 01=10=0,02=20=1,12=21=2 
    }
  }

  register int i;
  register site *s;
  char **local_pt[2][4];
  int mu, nu, index, gather, flip = 0, a, b;
  msg_tag *tag0[2], *tag1[2], *tag2[2], *tag3[2];

  for (mu = 0; mu < 4; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[4 + mu];
  }

  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(src[1], sizeof(matrix),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  tag1[0] = start_gather_site(F_OFFSET(link[0]), sizeof(matrix),
                              goffset[1], EVENANDODD, local_pt[0][1]);

  tag2[0] = start_gather_field(src[0], sizeof(matrix),
                               goffset[1], EVENANDODD, local_pt[0][2]);

  tag3[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                              goffset[0], EVENANDODD, local_pt[0][3]);

  // Main loop
  FORALLDIR(mu) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      index = plaq_index[mu][nu];
      gather = (flip + 1) % 2;
      if (index < NPLAQ - 1) {               // Start next set of gathers
        if (nu == NUMLINK - 1) {
          a = mu + 1;
          b = a + 1;
        }
        else {
          a = mu;
          b = nu + 1;
        }
        tag0[gather] = start_gather_field(src[b], sizeof(matrix), goffset[a],
                                          EVENANDODD, local_pt[gather][0]);

        tag1[gather] = start_gather_site(F_OFFSET(link[a]), sizeof(matrix),
                                         goffset[b], EVENANDODD,
                                         local_pt[gather][1]);

        tag2[gather] = start_gather_field(src[a], sizeof(matrix), goffset[b],
                                          EVENANDODD, local_pt[gather][2]);

        tag3[gather] = start_gather_site(F_OFFSET(link[b]), sizeof(matrix),
                                         goffset[a], EVENANDODD,
                                         local_pt[gather][3]);
      }

      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      wait_gather(tag2[flip]);
      wait_gather(tag3[flip]);
      FORALLSITES(i, s) {
        // Initialize dest[index][i]
        scalar_mult_nn(&(s->link[mu]), (matrix *)(local_pt[flip][0][i]),
                       s->bc[mu], &(plaq_dest[index][i]));

        // Add or subtract the other three terms
        mult_nn_dif(&(src[nu][i]), (matrix *)(local_pt[flip][1][i]),
                    &(plaq_dest[index][i]));

        scalar_mult_nn_dif(&(s->link[nu]), (matrix *)(local_pt[flip][2][i]),
                           s->bc[nu], &(plaq_dest[index][i]));

        mult_nn_sum(&(src[mu][i]), (matrix *)(local_pt[flip][3][i]),
                    &(plaq_dest[index][i]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      cleanup_gather(tag2[flip]);
      cleanup_gather(tag3[flip]);
      flip = gather;
    }
  }
}


void PtoLb(matrix *src[NPLAQ], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu;
  char **local_pt[2][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NPLAQ];

    for (mu = 0; mu < NPLAQ; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
  }
  
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};
  
  FORALLDIR(mu) {
  FORALLSITES(i, s) {           // edited ----- to initialize otherwise lead to NAN
    clear_mat(&(dest[mu][i])); } 
  } 
  
//--------  
 

  FORALLDIR(mu) {
  
        Real SIGN = SIGN_array[mu];

        tag0[mu] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                                         goffset[mu] , EVENANDODD,
                                         local_pt[0][mu]);             //phi(n+mu)
                                         
        tag1[mu] = start_gather_field(src[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // chi_mu(n+mu)                  

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
      
      
       mult_an(&(s->phi),(matrix *)(local_pt[1][mu][i]), &tmat );          // tmat = phi^bar(n)*chi_mu(n+mu)
       mult_na_dif( (matrix *)(local_pt[1][mu][i]),(matrix *)(local_pt[0][mu][i]), &tmat );      // tmat = phi^bar(n)*chi_mu(n+mu) - chi_mu(n+mu) phi(n+mu)
       scalar_mult_sum_matrix(&tmat,SIGN, &dest[mu][i]);    // dest[i] = dest[i] +SGN*{phi^bar(n)*chi_mu(n+mu) - chi_mu(n+mu) phi(n+mu)}
      
      

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
    
      
    }
   

 }


void LbtoP(matrix *src[NUMLINK], matrix *dest[NPLAQ]) {

  register int i;
  register site *s;
  int mu;
  char **local_pt[2][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NPLAQ];

    for (mu = 0; mu < NPLAQ; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
  }
  
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};
  
//--------  
 

  FORALLDIR(mu) {
  
        Real SIGN = SIGN_array[mu];

        tag0[mu] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[0][mu]);             //phi(n-mu)
                                         
        tag1[mu] = start_gather_field(src[mu], sizeof(matrix),
                                          goffset[mu] + 1, EVENANDODD,
                                          local_pt[1][mu]);              // psib_mu(n-mu)                  

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
      
      
       mult_an(&(s->phi),(matrix *)(local_pt[1][mu][i]), &tmat );          // tmat = phi^bar(n)*psib_mu(n-mu)
       mult_na_dif( (matrix *)(local_pt[1][mu][i]),(matrix *)(local_pt[0][mu][i]), &tmat );      // tmat = phi^bar(n)*psib_mu(n-mu) - psib_mu(n-mu) phi(n-mu)
       scalar_mult_sum_matrix(&tmat,SIGN, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{phi^bar(n)*psib_mu(n-mu) - psib_mu(n-mu) phi(n-mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
    
      
    }
   

 }

////////-------------


void PtoTb(matrix *src[NPLAQ], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu;
  char **local_pt[2][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NPLAQ];

    for (mu = 0; mu < NPLAQ; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
  }
  
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};
  
  FORALLDIR(mu) {
  FORALLSITES(i, s) {           // edited ----- to initialize otherwise lead to NAN
    clear_mat(&(dest[mu][i])); } 
  } 
  
//--------  
 

  FORALLDIR(mu) {
  
        Real SIGN = SIGN_array[mu];

        tag0[mu] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                                         goffset[mu] , EVENANDODD,
                                         local_pt[0][mu]);             //varphi(n+mu)
                                         
        tag1[mu] = start_gather_field(src[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // chi_mu(n+mu)                  

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
      
      
       mult_an(&(s->varphi),(matrix *)(local_pt[1][mu][i]), &tmat );          // tmat = varphi^bar(n)*chi_mu(n+mu)
       mult_na_dif( (matrix *)(local_pt[1][mu][i]),(matrix *)(local_pt[0][mu][i]), &tmat );      // tmat = varphi^bar(n)*chi_mu(n+mu) - chi_mu(n+mu) varphi(n+mu)
       scalar_mult_sum_matrix(&tmat,-1.0*SIGN, &dest[mu][i]);    // dest[i] = dest[i] +SGN*{varphi^bar(n)*chi_mu(n+mu) - chi_mu(n+mu) varphi(n+mu)}
      
      

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
    
      
    }
   

 }


void TbtoP(matrix *src[NUMLINK], matrix *dest[NPLAQ]) {

  register int i;
  register site *s;
  int mu;
  char **local_pt[2][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NPLAQ];

    for (mu = 0; mu < NPLAQ; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
  }
  
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};
  
//--------  
 

  FORALLDIR(mu) {
  
        Real SIGN = SIGN_array[mu];

        tag0[mu] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[0][mu]);             //varphi(n-mu)
                                         
        tag1[mu] = start_gather_field(src[mu], sizeof(matrix),
                                          goffset[mu] + 1, EVENANDODD,
                                          local_pt[1][mu]);              // psib_mu(n-mu)                  

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
      
      
       mult_an(&(s->varphi),(matrix *)(local_pt[1][mu][i]), &tmat );          // tmat = varphi^bar(n)*psib_mu(n-mu)
       mult_na_dif( (matrix *)(local_pt[1][mu][i]),(matrix *)(local_pt[0][mu][i]), &tmat );      // tmat = varphi^bar(n)*psib_mu(n-mu) - psib_mu(n-mu) varphi(n-mu)
       scalar_mult_sum_matrix(&tmat,-1.0*SIGN, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{varphi^bar(n)*psib_mu(n-mu) - psib_mu(n-mu) varphi(n-mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
    
      
    }
   

 }





#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat and tempmat2 for temporary storage
#ifdef VP
void Dminus(matrix *src[NPLAQ], matrix *dest[NUMLINK]) {

//edited 
      int p, q, indx;
 for (p = 0; p < NUMLINK; p++) {
    plaq_index[p][p] = -1;                                  // i,i=-1  ,  00=11=22=-1
    for (q = p + 1; q < NUMLINK; q++) {
      indx = p * (NUMLINK - 1) - p * (p + 1) / 2 + q - 1;//              01=0,02=1,12=2 
      plaq_index[p][q] = indx;
      plaq_index[q][p] = indx;                             // i,j=j,i, 01=10=0,02=20=1,12=21=2 
    }
  }

  register int i;
  register site *s;
  char **local_pt[2][2];
  int mu, nu, index, gather, flip = 0, a, b, next, opp_mu;
  matrix *mat[2];
  msg_tag *tag0[2], *tag1[2];

  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu];
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  // Start first set of gathers (mu = 1 and nu = 0)
  index = plaq_index[1][0];
  tag0[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                              goffset[0], EVENANDODD, local_pt[0][0]);

  FORALLSITES(i, s) {   // mu = 1 > nu = 0
    scalar_mult_nn(&(src[index][i]), &(s->link[1]), -1.0, &(mat[0][i]));
    FORALLDIR(mu)
      clear_mat(&(dest[mu][i]));        // Initialize
  }
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(nu) {
    FORALLDIR(mu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (nu < NUMLINK - 1 || mu < NUMLINK - 2) { // Start next set of gathers
        if (mu == NUMLINK - 1) {
          a = 0;
          b = nu + 1;
        }
        else if (mu == nu - 1) {
          a = mu + 2;
          b = nu;
        }
        else {
          a = mu + 1;
          b = nu;
        }
        next = plaq_index[a][b];
        tag0[gather] = start_gather_site(F_OFFSET(link[a]), sizeof(matrix),
                                         goffset[b], EVENANDODD,
                                         local_pt[gather][0]);

        FORALLSITES(i, s) {
          if (a > b) {      // src is anti-symmetric under a <--> b
            scalar_mult_nn(&(src[next][i]), &(s->link[a]), -1.0,
                           &(mat[gather][i]));
          }
          else {
            mult_nn(&(src[next][i]), &(s->link[a]), &(mat[gather][i]));
          }
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[a] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      opp_mu = OPP_LDIR(mu);
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu)      // src is anti-symmetric under mu <--> nu
          mult_nn_dif((matrix *)(local_pt[flip][0][i]), &(src[index][i]),
                      &(dest[nu][i]));
        else
          mult_nn_sum((matrix *)(local_pt[flip][0][i]), &(src[index][i]),
                      &(dest[nu][i]));

        scalar_mult_dif_matrix((matrix *)(local_pt[flip][1][i]),
                               s->bc[opp_mu], &(dest[nu][i]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
}
#endif
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
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(matrix), goffset[mu + 1],
                                       EVENANDODD, gen_pt[mu + 1]);

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


//zetab,psib dplus

void DplusSzbtoLb(matrix *src, matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;
 /* 
  FORALLDIR(mu) {
  FORALLSITES(i, s) {           // edited ----- to initialize otherwise lead to NAN
    clear_mat(&(dest[mu][i])); } 
  }  
*/
  tag[0] = start_gather_field(src, sizeof(matrix), goffset[0],  // gen_pt[0]=zb(i+0)
                              EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(matrix), goffset[mu + 1],
                                       EVENANDODD, gen_pt[mu + 1]);               // gen_pt[mu]=zb(i+mu)

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_nn( &(s->link[mu]),(matrix *)(gen_pt[mu][i]), &tmat);           // tmat = U_mu(i)zb(i+mu)
      //scalar_mult_matrix(&tmat, s->bc[mu], &tmat);                        
      mult_nn_dif(&(src[i]), &(s->link[mu]), &tmat);                       // tmat = U_mu(i)zb(i+mu) - zb(i)U_mu(i)
      scalar_mult_sum_matrix(&tmat, 0.5, &(dest[mu][i]));                  // dest[mu][i] += 0.5 { U_mu(i)zb(i+mu) - zn(i)U_mu(i) }
    }
    cleanup_gather(tag[mu]);
  }
}


//etab,thetab dplus

void DplusSebtoTb(matrix *src, matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;
  /*
  FORALLDIR(mu) {
  FORALLSITES(i, s) {           // edited ----- to initialize otherwise lead to NAN
    clear_mat(&(dest[mu][i])); } 
  }*/

  tag[0] = start_gather_field(src, sizeof(matrix), goffset[0],  // gen_pt[0]=zb(i+0)
                              EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(src, sizeof(matrix), goffset[mu + 1],
                                       EVENANDODD, gen_pt[mu + 1]);               // gen_pt[mu]=zb(i+mu)

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_nn( &(s->link[mu]),(matrix *)(gen_pt[mu][i]), &tmat);           // tmat = U_mu(i)zb(i+mu)
      //scalar_mult_matrix(&tmat, s->bc[mu], &tmat);                        
      mult_nn_dif(&(src[i]), &(s->link[mu]), &tmat);                       // tmat = U_mu(i)zb(i+mu) - zb(i)U_mu(i)
      scalar_mult_sum_matrix(&tmat, 0.5, &(dest[mu][i]));                  // dest[mu][i] += 0.5 { U_mu(i)zb(i+mu) - zn(i)U_mu(i) }
    }
    cleanup_gather(tag[mu]);
  }
}


#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from site source to link destination
//   U^{-1}[a](x) * sum_b {D[b][a](x) + D[a][b](x - b)}
// D is Tr[eta] * plaqdet, Tr[eta] = i sqrt(N) eta^D
// T is Tr[U^{-1} Lambda]
// Assume compute_plaqdet() has already been run
// Use tr_dest and tempdet for temporary storage
// bc[b](x - b) = bc[-b](x) on eta(x - b) psi_a(x)
// Add negative to dest instead of overwriting
// Negative sign is due to anti-commuting eta past psi
#ifdef SV
void detStoL(matrix *dest[NUMLINK]) {
  register int i;
  register site *s;
  int a, b, opp_b, next;
  Real localG = -0.5 * C2 * G;
  complex tc;
  msg_tag *tag[NUMLINK];

  // Save Tr[eta(x)] plaqdet[a][b](x)
  //   or Tr[eta(x)] ZWstar[a][b](x) in tempdet[a][b]
  FORALLDIR(a) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLSITES(i, s) {
        CMUL(tr_eta[i], plaqdet[a][b][i], tempdet[a][b][i]);
        CMUL(tr_eta[i], plaqdet[b][a][i], tempdet[b][a][i]);
      }
    }
  }

  // Now we gather tempdet in both cases
  // Start first gather for (a, b) = (0, 1)
  tag[1] = start_gather_field(tempdet[0][1], sizeof(complex),
                              goffset[1] + 1, EVENANDODD, gen_pt[1]);

  FORALLDIR(a) {
    // Initialize accumulator for sum over b
    FORALLSITES(i, s)
      tr_dest[i] = cmplx(0.0, 0.0);

    FORALLDIR(b) {
      if (a == b)
        continue;

      // Start next gather unless we're doing the last (a=4, b=3)
      next = b + 1;
      if (next < NUMLINK && a + b < 2 * NUMLINK - 3) {
        if (next == a)              // Next gather is actually (a, b + 2)
          next++;

        tag[next] = start_gather_field(tempdet[a][next], sizeof(complex),
                                       goffset[next] + 1, EVENANDODD,
                                       gen_pt[next]);
      }
      else if (next == NUMLINK) {   // Start next gather (a + 1, 0)
        tag[0] = start_gather_field(tempdet[a + 1][0], sizeof(complex),
                                    goffset[0] + 1, EVENANDODD, gen_pt[0]);
      }

      // Accumulate tempdet[b][a](x) + tempdet[a][b](x - b)
      opp_b = OPP_LDIR(b);
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        tr_dest[i].real += s->bc[opp_b] * tc.real;
        tr_dest[i].imag += s->bc[opp_b] * tc.imag;
        CSUM(tr_dest[i], tempdet[b][a][i]);
      }
      cleanup_gather(tag[b]);
    }

    // Multiply U_a^{-1} by sum, add to dest[a][i]
    FORALLSITES(i, s) {
      CMULREAL(tr_dest[i], localG, tc);
      c_scalar_mult_sum_mat(&(Uinv[a][i]), &tc, &(dest[a][i]));
    }
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
  register int i, mu, nu, opp_mu;
  register site *s;
  int gather = 1, flip = 0;
  matrix *mat[2];
  msg_tag *tag[NUMLINK];

  mat[0] = tempmat;
  mat[1] = tempmat2;

  FORALLSITES(i, s) {           // Set up first gather
    clear_mat(&(dest[i]));      // Initialize
    mult_an(&(s->link[0]), &(src[0][i]), &(mat[0][i])); //mat[0][i]=U_0^dag psi_0
  }
  tag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[0]=U_0^dag(i-0) psi_0(i-0)

  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_an(&(s->link[nu]), &(src[nu][i]), &(mat[gather][i]));  // mat[gather][i]= U_nu^dag(i)psi_nu(i)
      tag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); // gen_pt[nu]= U_nu^dag(i-nu)psi_nu(i-nu)
    }

    opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]), s->bc[opp_mu],
                             &(dest[i]));                                 // dest[i]= dest[i]-U_nu^dag(i-nu)psi_nu(i-nu)*bc[opp_mu]
      mult_na_sum(&(src[mu][i]), &(s->link[mu]), &(dest[i]));             // dest[i]= -U_nu^dag(i-nu)psi_nu(i-nu)*bc[opp_mu] + psi_mu(i)U_mu(i)^dag
    }
    cleanup_gather(tag[mu]);
    flip = gather;
  }

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix(&(dest[i]), 0.5, &(dest[i]));
}


// psib,zetabar Dminus

void DminusLbtoSzb(matrix *src[NUMLINK], matrix *dest) {
  register int i, mu, nu, opp_mu;
  register site *s;
  int gather = 1, flip = 0;
  matrix *mat[2];
  msg_tag *tag[NUMLINK];

  mat[0] = tempmat;
  mat[1] = tempmat2;

  FORALLSITES(i, s) {           // Set up first gather
    clear_mat(&(dest[i]));      // Initialize
    mult_nn(&(src[0][i]),&(s->link[0]), &(mat[0][i])); //mat[0][i]=psi_0 U_0
  }
  tag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[0]=psi_0(i-0)U_0(i-0) 

  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_nn(&(src[nu][i]),&(s->link[nu]), &(mat[gather][i]));  // mat[gather][i]= psi_nu(i)U_nu(i)
      tag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); // gen_pt[nu]= psi_nu(i-nu)U_nu(i-nu)
    }

   // opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]), 1.0,
                             &(dest[i]));                                 // dest[i]= dest[i]-psi_nu(i-nu)U_nu(i-nu)
      mult_nn_sum( &(s->link[mu]),&(src[mu][i]), &(dest[i]));             // dest[i]= dest[i]-psi_nu(i-nu)U_nu(i-nu) + U_mu(i)psi_mu(i)
    }
    cleanup_gather(tag[mu]);
    flip = gather;
  }

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix(&(dest[i]), 0.5, &(dest[i]));
}


// thetab, etab Dminus

void DminusTbtoSeb(matrix *src[NUMLINK], matrix *dest) {
  register int i, mu, nu, opp_mu;
  register site *s;
  int gather = 1, flip = 0;
  matrix *mat[2];
  msg_tag *tag[NUMLINK];

  mat[0] = tempmat;
  mat[1] = tempmat2;

  FORALLSITES(i, s) {           // Set up first gather
    clear_mat(&(dest[i]));      // Initialize
    mult_nn(&(src[0][i]),&(s->link[0]), &(mat[0][i])); //mat[0][i]=psi_0 U_0
  }
  tag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[0]=psi_0(i-0)U_0(i-0) 

  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_nn(&(src[nu][i]),&(s->link[nu]), &(mat[gather][i]));  // mat[gather][i]= psi_nu(i)U_nu(i)
      tag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); // gen_pt[nu]= psi_nu(i-nu)U_nu(i-nu)
    }

   // opp_mu = OPP_LDIR(mu);
    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]), 1.0,
                             &(dest[i]));                                 // dest[i]= dest[i]-psi_nu(i-nu)U_nu(i-nu)
      mult_nn_sum( &(s->link[mu]),&(src[mu][i]), &(dest[i]));             // dest[i]= dest[i]-psi_nu(i-nu)U_nu(i-nu) + U_mu(i)psi_mu(i)
    }
    cleanup_gather(tag[mu]);
    flip = gather;
  }

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix(&(dest[i]), 0.5, &(dest[i]));
}

#endif


#ifdef SV

void DplusSebtoSz(matrix *src, matrix *dest) {
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[0];
  matrix tmat;
  
  
  FORALLSITES(i, s) {           // edited ----- to initialize otherwise lead to NAN
    clear_mat(&(dest[i])); } 
    

  tag[0] = start_gather_field(src, sizeof(matrix), goffset[0],  // gen_pt[0]=eb(i+0)
                              EVENANDODD, gen_pt[0]);
 

    wait_gather(tag[0]);
    FORALLSITES(i, s) {
      mult_nn( &(s->varphi),(matrix *)(gen_pt[0][i]), &tmat);           // tmat = varphi(i)eb(i+0)
      //scalar_mult_matrix(&tmat, s->bc[mu], &tmat);                        
      mult_nn_dif(&(src[i]), &(s->varphi), &tmat);                       // tmat = varphi(i)eb(i+0) - eb(i)varphi(i)
      scalar_mult_sum_matrix(&tmat, 0.5, &(dest[i]));                  // dest[mu][i] += 0.5 { varphi(i)eb(i+0) - eb(i)varphi(i) }
    }
    cleanup_gather(tag[0]);
  
}

void DminusSztoSeb(matrix *src, matrix *dest) {
  register int i, mu, nu, opp_mu;
  register site *s;
  int gather = 1, flip = 0;
  matrix *mat[2];
  msg_tag *tag[0];

  mat[0] = tempmat;
  mat[1] = tempmat2;

  FORALLSITES(i, s) {           // Set up first gather
    clear_mat(&(dest[i]));      // Initialize
    mult_nn(&(src[i]),&(s->varphi), &(mat[0][i])); //mat[0][i]=z*varphi
  }
  tag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[0]=z(i-0)*varphi(i-0)



   // opp_mu = OPP_LDIR(mu);
    wait_gather(tag[0]);
    FORALLSITES(i, s) {
      scalar_mult_dif_matrix((matrix *)(gen_pt[0][i]), 1.0,
                             &(dest[i]));                                 // dest[i]= dest[i]-z(i-0)*varphi(i-0)
      mult_nn_sum( &(s->varphi),&(src[i]), &(dest[i]));             // dest[i]= dest[i]-z(i-0)*varphi(i-0) + varphi(i)z(i)
    }
    cleanup_gather(tag[0]);
    
  

  // Overall factor of 1/2
  FORALLSITES(i, s)
    scalar_mult_matrix(&(dest[i]), 0.5, &(dest[i]));
}





#endif


// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant coupling from link source to site destination
//   sum_{a, b} D[a][b](x) * {T[b](x) +  T[a](x + b)}
// D is plaqdet and T is Tr[U^{-1} psi]
// Assume compute_plaqdet() has already been run
// bc[b](x) on eta(x) psi_a(x + b)
// Use Tr_Uinv and tr_dest for temporary storage
// Add to dest instead of overwriting
// Has same sign as DbminusLtoS (negative comes from generator normalization)
#ifdef SV
void detLtoS(matrix *src[NUMLINK], matrix *dest) {
  register int i;
  register site *s;
  int a, b;
  Real localG = 0.5 * C2 * G * sqrt((Real)NCOL);
  complex tc, tc2;
  msg_tag *tag[NUMLINK];

  // Prepare Tr[U_a^{-1} psi_a] = sum_j Tr[U_a^{-1} Lambda^j] psi_a^j
  // and save in Tr_Uinv[a]
  FORALLSITES(i, s) {
    tr_dest[i] = cmplx(0.0, 0.0);   // Initialize
    FORALLDIR(a)
      Tr_Uinv[a][i] = complextrace_nn(&(Uinv[a][i]), &(src[a][i]));
  }

  // Start first gather of Tr[U_a^{-1} psi_a] from x + b for (0, 1)
  tag[1] = start_gather_field(Tr_Uinv[0], sizeof(complex),
                              goffset[1], EVENANDODD, gen_pt[1]);

  // Main loop
  FORALLDIR(a) {
    FORALLDIR(b) {
      if (a == b)
        continue;

      if (a == 0) {         // Start other gather (1, 0)
        tag[0] = start_gather_field(Tr_Uinv[1], sizeof(complex),
                                    goffset[0], EVENANDODD, gen_pt[0]);
      }

      // Accumulate D[a][b](x) {T[b](x) + T[a](x + b)} in tr_dest
      wait_gather(tag[b]);
      FORALLSITES(i, s) {
        tc = *((complex *)(gen_pt[b][i]));
        tc2.real = Tr_Uinv[b][i].real + s->bc[b] * tc.real;
        tc2.imag = Tr_Uinv[b][i].imag + s->bc[b] * tc.imag;
        CMUL(plaqdet[a][b][i], tc2, tc);
        // localG is purely imaginary...
        tr_dest[i].real -= tc.imag * localG;
        tr_dest[i].imag += tc.real * localG;
      }
      cleanup_gather(tag[b]);
    }
  }

  // Add to dest (negative comes from generator normalization)
  FORALLSITES(i, s)
    c_scalar_mult_dif_mat(&(Lambda[DIMF - 1]), &(tr_dest[i]), &(dest[i]));
}
#endif







#ifdef SS


//(1)

// Term1
// zeta-etab
// src = etab , dest = zeta   

void SebtoSz(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
  
   FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      } 

     FORALLSITES(i, s) {

       mult_nn(&(s->varphi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_nn_dif(&(src[i]), &(s->varphi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi(n)
       scalar_mult_sum_matrix(&tmat,1.0, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      }

}

//#endif

// Term2
// etab-zeta
// src = zeta , dest = etab

//#ifdef SS

void SztoSeb(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
/*  FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }    */                        // Initialize
            
     FORALLSITES(i, s) {

       mult_nn(&(s->varphi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_nn_dif(&(src[i]), &(s->varphi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi
       scalar_mult_sum_matrix(&tmat,1.0, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      } 

}


//(2)

// Term1
// zeta-ztab
// src = ztab , dest = zeta   

void SzbtoSz(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
  
 /*   FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }  */

     FORALLSITES(i, s) {

       mult_nn(&(s->phi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_nn_dif(&(src[i]), &(s->phi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi(n)
       scalar_mult_sum_matrix(&tmat,-1.0, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      }

}

//#endif

// Term2
// zetab-zeta
// src = zeta , dest = zetab

//#ifdef SS

void SztoSzb(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
/*  FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }     */                       // Initialize
            
     FORALLSITES(i, s) {

       mult_nn(&(s->phi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_nn_dif(&(src[i]), &(s->phi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi
       scalar_mult_sum_matrix(&tmat,-1.0, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      } 

}



//(3)

// Term1
// etab-eta
// src = etab , dest = eta   

void SebtoSe(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
  
/*   FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }  */

     FORALLSITES(i, s) {

       mult_an(&(s->phi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_na_dif(&(src[i]), &(s->phi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi(n)
       scalar_mult_sum_matrix(&tmat,0.5, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      }

}

//#endif

// Term2
// eta-etab
// src = eta , dest = etab

//#ifdef SS

void SetoSeb(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
 /* FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }     */                       // Initialize
            
     FORALLSITES(i, s) {

       mult_an(&(s->phi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_na_dif(&(src[i]), &(s->phi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi
       scalar_mult_sum_matrix(&tmat,0.5, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      } 

}


//(4)

// Term1
// zetab-eta
// src = zetab , dest = eta   

void SzbtoSe(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
  
 /*  FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }  */

     FORALLSITES(i, s) {

       mult_an(&(s->varphi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_na_dif(&(src[i]), &(s->varphi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi(n)
       scalar_mult_sum_matrix(&tmat,0.5, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      }

}

//#endif

// Term2
// eta-zetab
// src = eta , dest = zetab

//#ifdef SS

void SetoSzb(matrix *src, matrix *dest){
  register int i;
  register site *s;
  int mu;

  matrix tmat;
 /* FORALLSITES(i, s){
     clear_mat(&(dest[i])); 
      }     */                       // Initialize
            
     FORALLSITES(i, s) {

       mult_an(&(s->varphi), &(src[i]), &tmat );          // tmat = varphi(n)*etab(n)
       mult_na_dif(&(src[i]), &(s->varphi), &tmat );      // tmat = varphi(n)*etab(n) - etab(n)varphi
       scalar_mult_sum_matrix(&tmat,0.5, &dest[i]);    // dest[i] = dest[i] + 1.0*{ varphi(n)*etab(n) - etab(n)varphi }
       
      } 

}

#endif



#ifdef VV

//(1)

// Term1
// psi-psib
// src = psi[NUMLINK] , dest = psib[NUMLINK]   

void LtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]){
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;
  
/*  FORALLDIR(mu) {
   FORALLSITES(i, s){
     clear_mat(&(dest[mu][i])); 
      } 
    } */
    
    
  tag[0] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);  // gen_pt[0]= varphi(n+0)
    
  FORALLDIR(mu) {
     if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                              goffset[mu+1], EVENANDODD, gen_pt[mu+1]);  // gen_pt[mu]= varphi(n+mu)

    wait_gather(tag[mu]);
    
     FORALLSITES(i, s) {

       mult_nn(&(s->varphi), &(src[mu][i]), &tmat );          // tmat = varphi(n)*psi_mu(n)
       mult_nn_dif(&(src[mu][i]), (matrix *)(gen_pt[mu][i]), &tmat );      // tmat = varphi(n)*psi_mu(n) - psi_mu(n)varphi(n+mu)
       scalar_mult_sum_matrix(&tmat,-1.0, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{varphi(n)*psi_mu(n) - psi_mu(n)varphi(n+mu) }
       
      }
    cleanup_gather(tag[mu]);   
      
   }
}

//#endif

// Term2
// psib-psi
// src = psib[NUMLINK] , dest = psi[NUMLINK]   

void LbtoL(matrix *src[NUMLINK], matrix *dest[NUMLINK]){
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;
  
/*  FORALLDIR(mu) {
   FORALLSITES(i, s){
     clear_mat(&(dest[mu][i])); 
      } 
    }*/
    
    
  tag[0] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);  // gen_pt[0]= varphi(n+0)
    
  FORALLDIR(mu) {
     if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                              goffset[mu+1], EVENANDODD, gen_pt[mu+1]);  // gen_pt[mu]= varphi(n+mu)

    wait_gather(tag[mu]);
    
     FORALLSITES(i, s) {

       mult_nn((matrix *)(gen_pt[mu][i]), &(src[mu][i]), &tmat );          // tmat = varphi(n+i)*psib_mu(n)
       mult_nn_dif(&(src[mu][i]), &(s->varphi), &tmat );      // tmat = varphi(n+i)*psib_mu(n) - psib_mu(n)varphi(n)
       scalar_mult_sum_matrix(&tmat,-1.0, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{varphi(n+i)*psib_mu(n) - psib_mu(n)varphi(n)}
       
      }
    cleanup_gather(tag[mu]);   
      
   }
}

//(2) thetabar,psi


// Term1
// psi-thetab
// src = psi[NUMLINK] , dest = thetab[NUMLINK]   

void LtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]){
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;
  
/*  FORALLDIR(mu) {
   FORALLSITES(i, s){
     clear_mat(&(dest[mu][i])); 
      } 
    } */
    
    
  tag[0] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);  // gen_pt[0]= varphi(n+0)
    
  FORALLDIR(mu) {
     if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                              goffset[mu+1], EVENANDODD, gen_pt[mu+1]);  // gen_pt[mu]= varphi(n+mu)

    wait_gather(tag[mu]);
    
     FORALLSITES(i, s) {

       mult_nn(&(s->phi), &(src[mu][i]), &tmat );          // tmat = varphi(n)*psi_mu(n)
       mult_nn_dif(&(src[mu][i]), (matrix *)(gen_pt[mu][i]), &tmat );      // tmat = varphi(n)*psi_mu(n) - psi_mu(n)varphi(n+mu)
       scalar_mult_sum_matrix(&tmat,-1.0, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{varphi(n)*psi_mu(n) - psi_mu(n)varphi(n+mu) }
       
      }
    cleanup_gather(tag[mu]);   
      
   }
}

//#endif

// Term2
// thetab-psi
// src = thetab[NUMLINK] , dest = psi[NUMLINK]   

void TbtoL(matrix *src[NUMLINK], matrix *dest[NUMLINK]){
  register int i;
  register site *s;
  int mu;
  msg_tag *tag[NUMLINK];
  matrix tmat;
  
/*  FORALLDIR(mu) {
   FORALLSITES(i, s){
     clear_mat(&(dest[mu][i])); 
      } 
    }*/
    
    
  tag[0] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);  // gen_pt[0]= varphi(n+0)
    
  FORALLDIR(mu) {
     if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                              goffset[mu+1], EVENANDODD, gen_pt[mu+1]);  // gen_pt[mu]= varphi(n+mu)

    wait_gather(tag[mu]);
    
     FORALLSITES(i, s) {

       mult_nn((matrix *)(gen_pt[mu][i]), &(src[mu][i]), &tmat );          // tmat = varphi(n+i)*psib_mu(n)
       mult_nn_dif(&(src[mu][i]), &(s->phi), &tmat );      // tmat = varphi(n+i)*psib_mu(n) - psib_mu(n)varphi(n)
       scalar_mult_sum_matrix(&tmat,-1.0, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{varphi(n+i)*psib_mu(n) - psib_mu(n)varphi(n)}
       
      }
    cleanup_gather(tag[mu]);   
      
   }
}



// **** maybe the substraction  instead of addition the same problem which was with bosonic action


//---********----


 void DbplusTbtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[6][4];
  matrix tmat;
  msg_tag *tag0[6],*tag1[6],*tag2[6],*tag3[6];

  for (int b = 0; b < 6; b++) {
    for (int a = 0; a < 4; a++) {
        local_pt[b][a] = gen_pt[b * 4 + a];
    }
  }

 
 
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {-1,1,-1, 1,-1, 1};
  Real SGN;
  
     mu  = mu_vals[0];
     nu  = nu_vals[0];
     rho = rho_vals[0];
     
 
        tag0[0] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][0]);              // thetab_rho(n-rho)      
                                          
        tag1[0] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][1]);              // thetab_rho(n+mu)      
                                  

        tag2[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[0][2]);                // U_nu(n-nu)
                                         
        tag3[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[0][3]);                // U_nu(n+mu)  
 
 

  for (int j = 0; j < 6; j++) {
  
    
     
     

              
     if(j<5)

     { 
        
     mu  = mu_vals[j+1];
     nu  = nu_vals[j+1];
     rho = rho_vals[j+1];
     
     
  
        tag0[j+1] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[j+1][0]);              // thetab_rho(n-rho)      
                                          
        tag1[j+1] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[j+1][1]);              // thetab_rho(n+mu)      
                                  

        tag2[j+1] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[j+1][2]);                // U_nu(n-nu)
                                         
        tag3[j+1] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[j+1][3]);                // U_nu(n+mu)     
      }
      
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];      
     SGN = sgn_val[j];                                           
                                                                          
      wait_gather(tag0[j]);
      wait_gather(tag1[j]);
      wait_gather(tag2[j]);
      wait_gather(tag3[j]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[j][0][i]),(matrix *)(local_pt[j][3][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[j][2][i]),(matrix *)(local_pt[j][1][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,1.0*SGN, &dest[mu][i]);                             // dest[i] = dest[i] - 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag0[j]);
      cleanup_gather(tag1[j]);
      cleanup_gather(tag2[j]);
      cleanup_gather(tag3[j]);
    
     
              
      
    }
   

 }
 
 
 
 /////////////
 void DbminusLbtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[6][4];
  matrix tmat;
  msg_tag *tag0[6],*tag1[6],*tag2[6],*tag3[6];
 

  for (int b = 0; b < 6; b++) {
    for (int a = 0; a < 4; a++) {
        local_pt[b][a] = gen_pt[b * 4 + a];
    }
  }

 
 
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {1,-1, 1,-1, 1,-1};
  Real SGN;


  
     mu  = mu_vals[0];
     nu  = nu_vals[0];
     rho = rho_vals[0];
   

          
 
         
        tag0[0] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][0]);              // Lb_rho(n-rho)      
                                          
        tag1[0] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][1]);              // Lb_rho(n+mu)      
                                  

        tag2[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[0][2]);                // U_nu(n-nu)
                                         
        tag3[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[0][3]);                // U_nu(n+mu) 
  
 

  for (int j = 0; j < 6; j++) {
  
  
   if(j<5) 
    
   {  
  
     mu  = mu_vals[j+1];
     nu  = nu_vals[j+1];
     rho = rho_vals[j+1];
     
        

          
 
         
        tag0[j+1] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[j+1][0]);              // Lb_rho(n-rho)      
                                          
        tag1[j+1] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[j+1][1]);              // Lb_rho(n+mu)      
                                  

        tag2[j+1] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[j+1][2]);                // U_nu(n-nu)
                                         
        tag3[j+1] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[j+1][3]);                // U_nu(n+mu)     
                                         
    } 
    
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];    
     SGN = sgn_val[j];                                            
                                                                          
      wait_gather(tag0[j]);
      wait_gather(tag1[j]);
      wait_gather(tag2[j]);
      wait_gather(tag3[j]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[j][0][i]),(matrix *)(local_pt[j][3][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[j][2][i]),(matrix *)(local_pt[j][1][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,1.0*SGN, &dest[mu][i]);                        // dest[i] = dest[i] + 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag0[j]);
      cleanup_gather(tag1[j]);
      cleanup_gather(tag2[j]);
      cleanup_gather(tag3[j]);
      

      
 
  
            }
                  

 }




/*
 void DbplusTbtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho,a,gather,flip=0;
  char **local_pt[2][4];
  matrix tmat;
  msg_tag *tag0[2],*tag1[2],*tag2[2],*tag3[2];

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }
 
 
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {-1,1,-1, 1,-1, 1};
  Real SGN;
 

  for (int j = 0; j < 6; j++) {
  
     gather = (flip + 1) % 2;
     
     
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];
     SGN = sgn_val[j];
              


  
        tag0[flip] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[flip][0]);              // thetab_rho(n-rho)      
                                          
        tag1[flip] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[flip][1]);              // thetab_rho(n+mu)      
                                  

        tag2[flip] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[flip][2]);                // U_nu(n-nu)
                                         
        tag3[flip] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[flip][3]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      wait_gather(tag2[flip]);
      wait_gather(tag3[flip]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[flip][0][i]),(matrix *)(local_pt[flip][3][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[flip][2][i]),(matrix *)(local_pt[flip][1][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,1.0*SGN, &dest[mu][i]);                             // dest[i] = dest[i] - 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      cleanup_gather(tag2[flip]);
      cleanup_gather(tag3[flip]);
    
     // flip = gather;
              
      
    }
   

 }
 
 
 
 /////////////
 void DbminusLbtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho,a,gather,flip=0;
  char **local_pt[2][4];
  matrix tmat;
 // msg_tag *tag0[6],*tag1[6],*tag2[6],*tag3[6];
  msg_tag *tag0[2],*tag1[2],*tag2[2],*tag3[2];

  for (a = 0; a < 4; a++) {
    local_pt[0][a] = gen_pt[a];
    local_pt[1][a] = gen_pt[4 + a];
  }
 
 
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {1,-1, 1,-1, 1,-1};
  Real SGN;
  
 

  for (int j = 0; j < 6; j++) {
  
  
     gather = (flip + 1) % 2;
  
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];
     SGN = sgn_val[j];
        

          
 
         
        tag0[flip] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[flip][0]);              // Lb_rho(n-rho)      
                                          
        tag1[flip] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[flip][1]);              // Lb_rho(n+mu)      
                                  

        tag2[flip] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[flip][2]);                // U_nu(n-nu)
                                         
        tag3[flip] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[flip][3]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      wait_gather(tag2[flip]);
      wait_gather(tag3[flip]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[flip][0][i]),(matrix *)(local_pt[flip][3][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[flip][2][i]),(matrix *)(local_pt[flip][1][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,1.0*SGN, &dest[mu][i]);                        // dest[i] = dest[i] + 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      cleanup_gather(tag2[flip]);
      cleanup_gather(tag3[flip]);
      
    //  flip = gather;
      
 
  
            }
                  

 }


*/


/*

//---********----

 void DbplusTbtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
 
 
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {-1,1,-1, 1,-1, 1};
  Real SGN;
 

  for (int j = 0; j < 6; j++) {
  
  
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];
     SGN = sgn_val[j];
              


  
        tag[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // thetab_rho(n-rho)      
                                          
        tag[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // thetab_rho(n+mu)      
                                  

        tag[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag[mu]);
      wait_gather(tag[mu]);
      wait_gather(tag[mu]);
      wait_gather(tag[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                             // dest[i] = dest[i] - 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag[mu]);
      cleanup_gather(tag[mu]);
      cleanup_gather(tag[mu]);
      cleanup_gather(tag[mu]);
    
              
      
    }
   

 }
 
 
 
 /////////////
 void DbminusLbtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
  
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {1,-1, 1,-1, 1,-1};
  Real SGN;
  
 

  for (int j = 0; j < 6; j++) {
  
  
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];
     SGN = sgn_val[j];
        

          
 
         
        tag[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // Lb_rho(n-rho)      
                                          
        tag[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // Lb_rho(n+mu)      
                                  

        tag[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag[mu]);
      wait_gather(tag[mu]);
      wait_gather(tag[mu]);
      wait_gather(tag[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                        // dest[i] = dest[i] + 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag[mu]);
      cleanup_gather(tag[mu]);
      cleanup_gather(tag[mu]);
      cleanup_gather(tag[mu]);
      

      
 
  
            }
                  

 }
 
*/ 

/*
//---********----

 void DbplusTbtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NUMLINK], *tag2[NUMLINK], *tag3[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
 
 
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {-1,1,-1, 1,-1, 1};
  Real SGN;
 

  for (int j = 0; j < 6; j++) {
  
  
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];
     SGN = sgn_val[j];
              


  
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // thetab_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // thetab_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                             // dest[i] = dest[i] - 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);
    
              
      
    }
   

 }
 
 
 
 /////////////
 void DbminusLbtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NUMLINK], *tag2[NUMLINK], *tag3[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
  
  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val[6] = {1,-1, 1,-1, 1,-1};
  Real SGN;
  
 

  for (int j = 0; j < 6; j++) {
  
  
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];
     SGN = sgn_val[j];
        

          
 
         
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // Lb_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // Lb_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                        // dest[i] = dest[i] + 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);
      

      
 
  
            }
                  

 }
 
*/

/*

//---********----

 void DbplusTbtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NUMLINK], *tag2[NUMLINK], *tag3[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
  
  Real SGN;
  
//--------  
 

  FORALLDIR(mu) {
  

            
    //case1
            
    nu = (mu + 1) % 3;
    rho = (mu + 2) % 3;

    //mu = 0, nu = 1, rho = 2
    //mu = 1, nu = 2, rho = 0
    //mu = 2, nu = 0, rho = 1
    
    SGN = -1.0;
              


  
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // thetab_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // thetab_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                             // dest[i] = dest[i] - 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);
    

    // case 2
    nu = (mu + 2) % 3;
    rho = (mu + 1) % 3;

    //mu = 0, nu = 2, rho = 1
    //mu = 1, nu = 0, rho = 2
    //mu = 2, nu = 1, rho = 0      
           
    SGN = 1;       

          

  
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // thetab_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // thetab_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                            // dest[i] = dest[i] + 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);      
      

  
                       
                  
                  
      
    }
   

 }
 
 
 
 /////////////
 void DbminusLbtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NUMLINK], *tag2[NUMLINK], *tag3[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
  
  Real SGN;
  
//--------  
 

  FORALLDIR(mu) {
  
    //case1
            
    nu = (mu + 1) % 3;
    rho = (mu + 2) % 3;

    //mu = 0, nu = 1, rho = 2
    //mu = 1, nu = 2, rho = 0
    //mu = 2, nu = 0, rho = 1
    
    SGN = 1.0;
        

          
 
         
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // Lb_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // Lb_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                        // dest[i] = dest[i] + 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);
      
    // case 2
    nu = (mu + 2) % 3;
    rho = (mu + 1) % 3;

    //mu = 0, nu = 2, rho = 1
    //mu = 1, nu = 0, rho = 2
    //mu = 2, nu = 1, rho = 0      
           
    SGN = -1;       
      

         
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // Lb_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // Lb_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                        // dest[i] = dest[i] - 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]); 
      
 
  
            }
                  

 }
 
 
 
 */
 

/*

 void DbplusTbtoLb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NUMLINK], *tag2[NUMLINK], *tag3[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
  
  Real SGN;
  
//--------  
 

  FORALLDIR(mu) {
  
      FORALLDIR(nu) {
      
          if (nu == mu) continue;
      
      
          FORALLDIR(rho) {
          
          
          if (rho == mu || rho == nu) continue;
          
        
        // mu=i, nu=j,rho=k
        // mu=0,nu=1,rho=2
        // mu=0,nu=2,rho=1
        // mu=1,nu=0,rho=2
        // mu=1,nu=2,rho=0
        // mu=2,nu=0,rho=1
        // mu=2,nu=1,rho=0
          
        // Determine SGN using Levi-Civita parity
            if ((mu == 0 && nu == 1 && rho == 2) ||
                (mu == 1 && nu == 2 && rho == 0) ||
                (mu == 2 && nu == 0 && rho == 1)) {
                SGN = -1.0;
            } else {
                SGN = 1.0;
            }
  
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // thetab_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // thetab_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                                          // dest[i] = dest[i] - 1.0*{thetab_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)thetab_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);
  
                       }
                  
                  }
      
    }
   

 }


/////////////
 void DbminusLbtoTb(matrix *src[NUMLINK], matrix *dest[NUMLINK]) {

  register int i;
  register site *s;
  int mu,nu,rho;
  char **local_pt[4][3];
  matrix tmat;
  msg_tag *tag0[NUMLINK], *tag1[NUMLINK], *tag2[NUMLINK], *tag3[NUMLINK];

    for (mu = 0; mu < NUMLINK; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
    local_pt[2][mu] = gen_pt[6 + mu];
    local_pt[3][mu] = gen_pt[9 + mu];
  }
  
  Real SGN;
  
//--------  
 

  FORALLDIR(mu) {
  
      FORALLDIR(nu) {
      
          if (nu == mu) continue;
      
      
          FORALLDIR(rho) {
          
          
          if (rho == mu || rho == nu) continue;
          
        
        // mu=k, nu=j,rho=i 
        // mu=0,nu=1,rho=2
        // mu=0,nu=2,rho=1
        // mu=1,nu=0,rho=2
        // mu=1,nu=2,rho=0
        // mu=2,nu=0,rho=1
        // mu=2,nu=1,rho=0
          
        // Determine SGN using Levi-Civita parity
            if ((mu == 0 && nu == 1 && rho == 2) ||
                (mu == 1 && nu == 2 && rho == 0) ||
                (mu == 2 && nu == 0 && rho == 1)) {
                SGN = 1.0;                              // sign changes on flip theta <-> linkbar
            } else {
                SGN = -1.0;
            }
  
        tag0[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][mu]);              // Lb_rho(n-rho)      
                                          
        tag1[mu] = start_gather_field(src[rho], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // Lb_rho(n+mu)      
                                  

        tag2[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),  
                                         goffset[nu] + 1, EVENANDODD,
                                         local_pt[2][mu]);                // U_nu(n-nu)
                                         
        tag3[mu] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                         goffset[mu], EVENANDODD,
                                         local_pt[3][mu]);                // U_nu(n+mu)             
                                                                          
      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      wait_gather(tag2[mu]);
      wait_gather(tag3[mu]);
      
      FORALLSITES(i, s) {
       
       mult_na((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[3][mu][i]), &tmat );          // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu)
       mult_an_dif((matrix *)(local_pt[2][mu][i]),(matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu) 
       scalar_mult_sum_matrix(&tmat,SGN, &dest[mu][i]);                                        // dest[i] = dest[i] - 1.0*{Lb_rho(n-rho)*Ubar_nu(n+mu) - Ubar_nu(n-nu)Lb_rho(n+mu)}

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
      cleanup_gather(tag2[mu]);
      cleanup_gather(tag3[mu]);
  
                       }
                  
                  }
      
    }
   

 }


*/

#endif




#ifdef SP
void DbplusPtoSz(matrix *src[NPLAQ], matrix *dest) {

  register int i;
  register site *s;
  int mu;
  char **local_pt[2][3];
  matrix tmat;
  msg_tag *tag0[NPLAQ], *tag1[NPLAQ];

    for (mu = 0; mu < NPLAQ; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
  }
  
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};
  
//--------  
 

  FORALLDIR(mu) {
  
        Real SIGN = SIGN_array[mu];

        tag0[mu] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[0][mu]);             //U_a(n-a)
                                         
        tag1[mu] = start_gather_field(src[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);              // chi_mu(n+mu)                  

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
      
      
       mult_na((matrix *)(local_pt[1][mu][i]),&(s->link[mu]), &tmat );          // tmat = chi_mu(n+mu) *Ubar_mu(n)
       mult_an_dif( (matrix *)(local_pt[0][mu][i]), &(src[mu][i]), &tmat );      // tmat = chi_mu(n+mu) *Ubar_mu(n) - Ubar_mu(n-mu)chi_mu(n)
       scalar_mult_sum_matrix(&tmat,SIGN, &dest[i]);    // dest[i] = dest[i] - 1.0*{chi_mu(n+mu) *Ubar_mu(n) - Ubar_mu(n-mu)chi_mu(n)}
      
      

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
  
      
      
    }
   
 
 
  
  
  
 }
 
 
 void DbminusSztoP(matrix *src, matrix *dest[NPLAQ]) {

  register int i;
  register site *s;
  int mu;
  char **local_pt[2][3];
  matrix tmat;
  msg_tag *tag0[NPLAQ], *tag1[NPLAQ];

    for (mu = 0; mu < NPLAQ; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[3 + mu];
  }
  
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};
  
//--------  
 

  FORALLDIR(mu) {
  
        Real SIGN = SIGN_array[mu];

        tag0[mu] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[0][mu]);             //U_a(n-a)
                                         
        tag1[mu] = start_gather_field(src, sizeof(matrix),
                                          goffset[mu]+1, EVENANDODD,
                                          local_pt[1][mu]);              // zeta(n-mu)                  

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
      
      
       mult_na(&(src[i]),(matrix *)(local_pt[0][mu][i]), &tmat );          // tmat = zeta(n)*Ubar_mu(n-mu)
       mult_an_dif( (matrix *)(local_pt[0][mu][i]), (matrix *)(local_pt[1][mu][i]), &tmat );      // tmat = zeta(n)*Ubar_mu(n-mu) - Ubar_mu(n-mu)zeta(n-mu)
       scalar_mult_sum_matrix(&tmat,SIGN, &dest[mu][i]);    // dest[i] = dest[i] - 1.0*{zeta(n)*Ubar_mu(n-mu) - Ubar_mu(n-mu)zeta(n-mu)}
      
      

      }
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
  
      
      
    }
   
 
 
  
  
  
 }
 

#endif





// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Twist_Fermion matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
void fermion_op(Twist_Fermion *src, Twist_Fermion *dest, int sign) {
  register int i, mu;
  register site *s;

  // Copy src TwistFermion into fieldwise site, link and plaq fermions,
  // overwriting all of the latter
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(src[i].Fsite), &(site_src[i]));
      mat_copy(&(src[i].Fsitezb), &(sitezb_src[i]));
      mat_copy(&(src[i].Fsiteeb), &(siteeb_src[i]));
      mat_copy(&(src[i].Fsitez), &(sitez_src[i]));
      FORALLDIR(mu){
        mat_copy(&(src[i].Flink[mu]), &(link_src[mu][i]));
        mat_copy(&(src[i].Flinkb[mu]), &(linkb_src[mu][i]));
        mat_copy(&(src[i].Fthetab[mu]), &(thetab_src[mu][i]));
        }
      for (mu = 0; mu < NPLAQ; mu++)
        mat_copy(&(src[i].Fplaq[mu]), &(plaq_src[mu][i]));
    }
  }
  else if (sign == -1) {
    FORALLSITES(i, s) {
      adjoint(&(src[i].Fsite), &(site_src[i]));
      adjoint(&(src[i].Fsitezb), &(sitezb_src[i]));
      adjoint(&(src[i].Fsiteeb), &(siteeb_src[i]));
      adjoint(&(src[i].Fsitez), &(sitez_src[i]));
      FORALLDIR(mu){
        adjoint(&(src[i].Flink[mu]), &(link_src[mu][i]));
        adjoint(&(src[i].Flinkb[mu]), &(linkb_src[mu][i]));
        adjoint(&(src[i].Fthetab[mu]), &(thetab_src[mu][i]));
        }
      for (mu = 0; mu < NPLAQ; mu++)
        adjoint(&(src[i].Fplaq[mu]), &(plaq_src[mu][i]));
    }
  }
  else {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }
  FORALLSITES(i, s)
    tr_eta[i] = trace(&(site_src[i]));

  // Assemble separate routines for each term in the fermion operator
#ifdef VP
  Dplus(link_src, plaq_dest);             // Overwrites plaq_dest
  Dminus(plaq_src, link_dest);            // Overwrites link_dest
  
  PtoLb(plaq_src,linkb_dest);             // Overwrites linkb_dest 
  LbtoP(linkb_src,plaq_dest);             // Adds to plaq_dest
  
  PtoTb(plaq_src,thetab_dest);             // Overwrites thetab_dest 
  TbtoP(thetab_src,plaq_dest);             // Adds to plaq_dest
  
  
#else
  FORALLSITES(i, s) {                     // Zero link_dest and plaq_dest
    FORALLDIR(mu)
      clear_mat(&(link_dest[mu][i]));
    for (mu = 0; mu < NPLAQ; mu++)
      clear_mat(&(plaq_dest[mu][i]));
  }
#endif

#ifdef SV
  DbplusStoL(site_src, link_dest);        // Adds to link_dest

  // Site-to-link plaquette determinant contribution if G is non-zero
  // Only depends on Tr[eta(x)]
  if (doG)
    detStoL(link_dest);                   // Adds to link_dest

  DbminusLtoS(link_src, site_dest);       // Overwrites site_dest

  // Link-to-site plaquette determinant contribution if G is non-zero
  if (doG)
    detLtoS(link_src, site_dest);         // Adds to site_dest
    
  DplusSzbtoLb(sitezb_src, linkb_dest);        // overwrites -> Adds to linkb_dest
  DminusLbtoSzb(linkb_src, sitezb_dest);       // Overwrites sitezb_dest  
  
  DplusSebtoTb(siteeb_src, thetab_dest);        // overwrites -> Adds to thetab_dest
  DminusTbtoSeb(thetab_src, siteeb_dest);       // Overwrites siteeb_dest
  
  
#else
  FORALLSITES(i, s)                       // Zero site_dest
    clear_mat(&(site_dest[i]));
  FORALLSITES(i, s)                       // Zero site_dest
    clear_mat(&(sitezb_dest[i]));   
#endif


#ifdef SS


//1
  SebtoSz(siteeb_src, sitez_dest);        // Overwrites sitez_dest
  SztoSeb(sitez_src, siteeb_dest);        // Adds to siteeb_dest 
 
//2  
  SzbtoSz(sitezb_src, sitez_dest);        // Adds to sitez_dest
  SztoSzb(sitez_src, sitezb_dest);        // Add to sitezb_dest
  
//3  
  SebtoSe(siteeb_src, site_dest);        // Add to site_dest
  SetoSeb(site_src, siteeb_dest);        // Add to siteeb_dest 
//4  
  SzbtoSe(sitezb_src, site_dest);        // Add to site_dest
  SetoSzb(site_src, sitezb_dest);        // Add to sitezb_dest 
  
#else
/*
  FORALLSITES(i, s)   {    
  //  clear_mat(&(site_dest[i]));  
    clear_mat(&(siteeb_dest[i]));         // Zero siteeb_dest
    clear_mat(&(sitez_dest[i]));       
    clear_mat(&(sitezb_dest[i]));         // Zero sitezb_dest         
    
    
    } */
#endif 


#ifdef VV

LtoLb(link_src, linkb_dest);               // Adds to linkb_dest
LbtoL(linkb_src, link_dest);               // Adds to link_dest

LtoTb(link_src, thetab_dest);               // Adds to thetab_dest
TbtoL(thetab_src, link_dest);               // Adds to link_dest  

DbplusTbtoLb(thetab_src,linkb_dest);       // Adds to linkb_dest
DbminusLbtoTb(linkb_src,thetab_dest);      // Adds to thetab_dest

#else
 
 /*   
 FORALLSITES(i, s) {                     // Zero linkb_dest 
    FORALLDIR(mu)
      clear_mat(&(linkb_dest[mu][i]));
  }
     
 FORALLSITES(i, s) {                     // Zero link_dest 
    FORALLDIR(mu)
      clear_mat(&(thetab_dest[mu][i]));
  } 
    */
#endif 


#ifdef SP

DbplusPtoSz(plaq_src,sitez_dest);        // Adds to sitez_dest
DbminusSztoP(sitez_src,plaq_dest);       // Adds to plaq_dest

#else

#endif


  // Copy local plaquette, link and site fermions into dest TwistFermion
  if (sign == 1) {
    FORALLSITES(i, s) {
      mat_copy(&(site_dest[i]), &(dest[i].Fsite));
      mat_copy(&(sitezb_dest[i]), &(dest[i].Fsitezb));
      mat_copy(&(siteeb_dest[i]), &(dest[i].Fsiteeb));
      mat_copy(&(sitez_dest[i]), &(dest[i].Fsitez));
      FORALLDIR(mu){
        mat_copy(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
        mat_copy(&(linkb_dest[mu][i]), &(dest[i].Flinkb[mu]));
        mat_copy(&(thetab_dest[mu][i]), &(dest[i].Fthetab[mu]));
        }
      for (mu = 0; mu < NPLAQ; mu++)
        mat_copy(&(plaq_dest[mu][i]), &(dest[i].Fplaq[mu]));
    }
  }
  else if (sign == -1) {    // Both negate and conjugate
    FORALLSITES(i, s) {
      neg_adjoint(&(site_dest[i]), &(dest[i].Fsite));
      neg_adjoint(&(sitezb_dest[i]), &(dest[i].Fsitezb));
      neg_adjoint(&(siteeb_dest[i]), &(dest[i].Fsiteeb));
      neg_adjoint(&(sitez_dest[i]), &(dest[i].Fsitez));
      FORALLDIR(mu){
        neg_adjoint(&(link_dest[mu][i]), &(dest[i].Flink[mu]));
        neg_adjoint(&(linkb_dest[mu][i]), &(dest[i].Flinkb[mu]));
        neg_adjoint(&(thetab_dest[mu][i]), &(dest[i].Fthetab[mu]));
        }
      for (mu = 0; mu < NPLAQ; mu++)
        neg_adjoint(&(plaq_dest[mu][i]), &(dest[i].Fplaq[mu]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Squared Twist_Fermion matrix--vector operation
//   dest = (D^2 + fmass^2).src
// Use tempTF for temporary storage
void DSq(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;

  fermion_op(src, tempTF, PLUS);
  fermion_op(tempTF, dest, MINUS);
  if (fmass > IMAG_TOL) {           // Assume fmass non-negative
    Real fmass2 = fmass * fmass;
    FORALLSITES(i, s)
      scalar_mult_sum_TF(&(src[i]), fmass2, &(dest[i]));
  }
}
// -----------------------------------------------------------------
