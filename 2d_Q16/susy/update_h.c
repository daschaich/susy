// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
#include "susy_includes.h"
// -----------------------------------------------------------------

//-----edited--------

// Set up translation of (i, j) to linear index of anti-symmetric matrix


//void setup_plaq_index() {




//}


// -----------------------------------------------------------------
// Update mom[NUMLINK] with the gauge force
// Include tunable coefficient C2 in the d^2 term of the action
// Use tr_dest, tempmat and tempdet for temporary storage
// Assume compute_plaqdet(), compute_DmuUmu()
// and compute_Fmunu() have already been run

int count;

double gauge_force(Real eps) {

  count=count+1;
 // printf("No of time pass = %d\n",count);

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

  register int i, mu, nu;
  register site *s;
  char **local_pt[2][2];
  //---edited----- char **local_pt[2][2];
  //char **local_pt[2][3]; 
  int a, b, gather, flip = 0, index, next;
  double returnit = 0.0, tr;
  complex tc,tcl;
  matrix tmat, tmat2, *mat[2];
  msg_tag *tag[NUMLINK], *tag0[2], *tag1[2];
  //---edited----- *tag0[2], *tag1[2]
  
  // msg_tag *tag[NUMLINK], *tag0[3], *tag1[3]; 



  // All contributions from d^2 term need a factor of C2
  // First we have the finite difference operator derivative times DmuUmu
  // Ubar_a(x) DmuUmu(x) - DmuUmu(x + a) Ubar_a(x)
  tag[0] = start_gather_field(DmuUmu, sizeof(matrix),
                              goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1)     // Start next gather
      tag[mu + 1] = start_gather_field(DmuUmu, sizeof(matrix),
                                       goffset[mu + 1], EVENANDODD,
                                       gen_pt[mu + 1]);

    wait_gather(tag[mu]);
    FORALLSITES(i, s) {
      mult_an(&(s->link[mu]), &(DmuUmu[i]), &(s->f_U[mu]));   // Initialize , ->force overwritten
      mult_na_dif((matrix *)(gen_pt[mu][i]), &(s->link[mu]), &(s->f_U[mu]));
    }
    cleanup_gather(tag[mu]);
  // printing mu   
  //  printf("mu = %d\n",mu);
  }





  // Next we have the plaquette determinant derivative contribution
  //   U_mu^{-1}(x) 2G sum_nu {D[nu][mu](x) + D[mu][nu](x-nu)}
  // D is Tr[DmuUmu] plaqdet[mu][nu], saved in tempdet[mu][nu]
  // Only compute if G is non-zero
  // Use tr_dest for temporary storage
  if (doG) {
    FORALLSITES(i, s) {
      tc = trace(&DmuUmu[i]);
      FORALLDIR(mu) {
        for (nu = mu + 1; nu < NUMLINK; nu++) {
          CMUL(tc, plaqdet[mu][nu][i], tempdet[mu][nu][i]);
          CMUL(tc, plaqdet[nu][mu][i], tempdet[nu][mu][i]);
        }
      }
    }

    // Start first gather (mu = 0 and nu = 1), labelled by nu
    // Gather D[mu][nu] from x - nu
    tag[1] = start_gather_field(tempdet[0][1], sizeof(complex),
                                goffset[1] + 1, EVENANDODD, gen_pt[1]);

    // Main loop
    FORALLDIR(mu) {
      // Zero tr_dest to hold sum
      FORALLSITES(i, s)
        tr_dest[i] = cmplx(0.0, 0.0);

      FORALLDIR(nu) {
        if (mu == nu)
          continue;

        if (mu < NUMLINK - 1 || nu < NUMLINK - 2) { // Start next gather
          if (nu == NUMLINK - 1) {
            a = mu + 1;
            b = 0;
          }
          else if (nu == mu - 1) {
            a = mu;
            b = nu + 2;
          }
          else {
            a = mu;
            b = nu + 1;
          }
          tag[b] = start_gather_field(tempdet[a][b], sizeof(complex),
                                      goffset[b] + 1, EVENANDODD,
                                      gen_pt[b]);
        }

        // Add D[nu][mu](x) to sum while gather runs
        FORALLSITES(i, s)
          CSUM(tr_dest[i], tempdet[nu][mu][i]);

        // Add D[mu][nu](x - nu) to sum
        wait_gather(tag[nu]);
        FORALLSITES(i, s)
          CSUM(tr_dest[i], *((complex *)(gen_pt[nu][i])));
        cleanup_gather(tag[nu]);
      }

      // Now add to force
      FORALLSITES(i, s) {
        CMULREAL(tr_dest[i], G, tc);
        c_scalar_mult_sum_mat(&(Uinv[mu][i]), &tc, &(s->f_U[mu]));
      }
    }
  }

  // Overall factor of C2 on all d^2 contributions
  if (C2 - 1.0 > IMAG_TOL) {
    FORALLSITES(i, s) {
      FORALLDIR(mu)
        scalar_mult_matrix(&(s->f_U[mu]), C2, &(s->f_U[mu]));
    }
  }





  // Contribution from Fbar_{ab} F_{ab} term
  // Start first set of gathers (mu = 0 and nu = 1)
  
  //2[U_b(n + µ_c )Fbar_cb (n) + Fbar_ac (n − µ_a )U_a (n − µ_a)]^T
  

  
  for (mu = 0; mu < 2; mu++) {
    local_pt[0][mu] = gen_pt[mu];
    local_pt[1][mu] = gen_pt[2 + mu]; //
  }
  mat[0] = tempmat;
  mat[1] = tempmat2;

  tag0[0] = start_gather_site(F_OFFSET(link[1]), sizeof(matrix),
                              goffset[0], EVENANDODD, local_pt[0][0]); // local_pt[0][0]=U_1(n+mu_0)

  index = plaq_index[0][1]; // index = 0
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_an(&(Fmunu[index][i]), &(s->link[1]), &(mat[0][i])); // mat[0]=F^bar[0,1](n)U_1(n)
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]); // // local_pt[0][1]=F^bar[0,1](n-mu_1)U_1(n-mu_1)

  // Main loop
  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;
     
    // node0_printf("mu,nu = %d,%d \n",mu,nu);
     
      gather = (flip + 1) % 2; // gather = remainder 1,0,1...
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) { // Start next gathers , mu=0,1 or nu=0 => [mu,nu]=01,02,10,12, or 10,20
        if (nu == NUMLINK - 1) {  // munu=02,12
          a = mu + 1;             // a=1,2
          b = 0;                  //b=0 ,ab=10,20
        }
        else if (nu == mu - 1) { //10,
          a = mu;                //a=1
          b = nu + 2;            //b=2  ,ab=12
        }
        else {                   //munu=01,20
          a = mu;                //a=0,2
          b = nu + 1;            //b=2,1 ,ab=02,01,22,21
        }

     //  node0_printf("a,b = %d,%d \n",mu,nu);

        tag0[gather] = start_gather_site(F_OFFSET(link[b]),
                                         sizeof(matrix), goffset[a],
                                         EVENANDODD, local_pt[gather][0]); // local_pt[gather][0]=U_b(n+a)

        next = plaq_index[a][b];
        
     //   node0_printf("a,b = %d,%d and indexab = %d \n",a,b,index);
        
        FORALLSITES(i, s) {
          if (a > b) {
            scalar_mult_matrix(&(Fmunu[next][i]), -1.0, &tmat); // tmat=-F[index]
            mult_an(&tmat, &(s->link[b]), &(mat[gather][i]));  // mat[gather][i]=-Fbar[index](n)U_b(n)
          }
          else
            mult_an(&(Fmunu[next][i]), &(s->link[b]), &(mat[gather][i])); //  mat[gather][i]=Fbar[index](n)U_b(n)
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);         //  local_pt[gather][1]=Fbar[index](n-mu_b)U_b(n-mu_b)
      }
      
      

      index = plaq_index[mu][nu];
      
      
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {
          scalar_mult_matrix(&(Fmunu[index][i]), -1.0, &tmat);  // tmat=-F[index](n)
          mult_na((matrix *)local_pt[flip][0][i], &tmat, &tmat2); // tmat2=[local_pt[gather][0]=U_b(n+a)]X-Fbar[index](n)
        }
        else
          mult_na((matrix *)local_pt[flip][0][i], &(Fmunu[index][i]), &tmat2); // tmat2=[ local_pt[flip][0]=U_b(n+a) ] X +Fbar[index](n)

        sub_matrix(&tmat2, (matrix *)local_pt[flip][1][i], &tmat); // tmat = [ local_pt[flip][0]=U_b(n+a) ] X +Fbar[index](n) - [ local_pt[gather][1]=Fbar[index](n-mu_b)U_b(n-mu_b) ]
        scalar_mult_sum_matrix(&tmat, 2.0, &(s->f_U[mu])); // f_U[mu] = f_U[mu] + 2* ( [ local_pt[flip][0]=U_b(n+a) ] X +Fbar[index](n) - [ local_pt[gather][1]=Fbar[index](n-mu_b)U_b(n-mu_b) ] )
        

      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
 


 
 // Here to put edited.......................
 
 
 
 // (1.0) contribution from DmuUmu*phi_commut (DmuUmu,phi_commut defined in action.c) wrt link

 // --------------------------- edited 7/9/2023------------------------------
  
   
  FORALLDIR(mu) {
  
   tag0[0] = start_gather_field(phi_commut, sizeof(matrix), // start_gather_site/field ?????????????????? field -> pointer, phi_commute defined as pointer in lattice.h
   
                              goffset[mu], EVENANDODD, gen_pt[0]);// gen_pt[0] = phi_commut(n+mu)


    wait_gather(tag0[0]);
    FORALLSITES(i, s) {
    
      mult_an(&(s->link[mu]), &(phi_commut[i]), &tmat);   // Initialize c= adag*b => tmat_a[n] = U^dag_a(n)PC(n) // a = mu
      mult_na_dif((matrix *)(gen_pt[0][i]), &(s->link[mu]), &tmat); //c=c-a*bdag => tmat_a[n] = U^dag_a(n)PC(n) - PC(n+a)U^dag_a(n)
      scalar_mult_sum_matrix(&tmat, 1.0, &(s->f_U[mu])); // c = c + s * b => f_U[mu] = f_U[mu] + 1.0*U^dag_a(n)PC(n) - PC(n+a)U^dag_a(n)
      
   
  }
  
   cleanup_gather(tag0[0]);
  
  }
  
  

 
  // (2.0) contribution from DmuUmu*varphi_commut (DmuUmu,varphi_commut defined in action.c) wrt link
  
 FORALLDIR(mu) {
  
   tag0[0] = start_gather_field(varphi_commut, sizeof(matrix), // start_gather_site/field ?????????????????? field -> pointer, varphi_commute defined as pointer in lattice.h
   
                              goffset[mu], EVENANDODD, gen_pt[0]);// gen_pt[0] = varphi_commut(n+mu)


    wait_gather(tag0[0]);
    FORALLSITES(i, s) {
    
      mult_an(&(s->link[mu]), &(varphi_commut[i]), &tmat);   // Initialize c= adag*b => tmat_a[n] = U^dag_a(n)PC(n) // a = mu
      mult_na_dif((matrix *)(gen_pt[0][i]), &(s->link[mu]), &tmat); //c=c-a*bdag => tmat_a[n] = U^dag_a(n)PC(n) - PC(n+a)U^dag_a(n)
      scalar_mult_sum_matrix(&tmat, 1.0, &(s->f_U[mu])); // c = c + s * b => f_U[mu] = f_U[mu] + 1.0*U^dag_a(n)PC(n) - PC(n+a)U^dag_a(n)
      
    
   
  }
  
   cleanup_gather(tag0[0]);
  
  }
  

  
  

  // D+ original
    
  //(3.0) contribution from -2|DmuPhi|^2 term wrt link
                                                             
                                                                                    
  FORALLDIR(mu) {
  
   tag0[0] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                              goffset[mu], EVENANDODD, gen_pt[0]);// gen_pt[0] = phi(n+mu)


    wait_gather(tag0[0]);
    
    FORALLSITES(i, s) {
      mult_na((matrix *)(gen_pt[0][i]), &(Dmuphi[mu][i]), &tmat);   // Initialize c= adag*b => tmat = phi_(i+mu)Dmuphi^dag[mu](i) // a = mu
      mult_an_dif(&(Dmuphi[mu][i]), &(s->phi), &tmat); //c=c-adag*b => tmat = phi_(i+mu)Dmuphi^dag[mu](i) - Dmuphi[mu]^dag(i)phi(i)
      scalar_mult_sum_matrix(&tmat,2.0, &(s->f_U[mu])); // c = c + s * b => f_U[mu] = f_U[mu]+ phi_(i+mu)Dmuphi^dag[mu](i) - Dmuphi[mu]^dag(i)phi(i) 
     } 
     
      cleanup_gather(tag0[0]);
  } 
  
  
 

  //(4.0) contribution from -2|Dmuvarphi|^2 term wrt link
                                                             
                                                                                    
  FORALLDIR(mu) {
  
   tag0[0] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                              goffset[mu], EVENANDODD, gen_pt[0]);// gen_pt[0] = phi(n+mu)


    wait_gather(tag0[0]);
    
    FORALLSITES(i, s) {
      mult_na((matrix *)(gen_pt[0][i]), &(Dmuvarphi[mu][i]), &tmat);   // Initialize c= adag*b => tmat = phi_(i+mu)Dmuphi^dag[mu](i) // a = mu
      mult_an_dif(&(Dmuvarphi[mu][i]), &(s->varphi), &tmat); //c=c-adag*b => tmat = phi_(i+mu)Dmuphi^dag[mu](i) - Dmuphi[mu]^dag(i)phi(i)
      scalar_mult_sum_matrix(&tmat,2.0, &(s->f_U[mu])); // c = c + s * b => f_U[mu] = f_U[mu]+ phi_(i+mu)Dmuphi^dag[mu](i) - Dmuphi[mu]^dag(i)phi(i) 
     } 
     
      cleanup_gather(tag0[0]);
  } 
  
  

 
   

// Only compute susy-breaking scalar potential term if bmass non-zero
  if (bmass > IMAG_TOL) {
    Real dmu;
#ifdef EIG_POT
    dmu =  2.0 * bmass * bmass;
#else
    Real tr;
    dmu = 2.0 * one_ov_N * bmass * bmass;
   
   //dmu = 0.0; // edited
#endif

    FORALLSITES(i, s) {
      FORALLDIR(mu) {
#ifdef EIG_POT
        // Ubar_a(x) [U_a(x) Ubar_a(x) - I]
        
        mult_na(&(s->link[mu]), &(s->link[mu]), &tmat);
        scalar_add_diag(&tmat, -1.0);
        scalar_mult_an_sum(&(s->link[mu]), &tmat, dmu, &(s->f_U[mu]));
        
#else
        // Ubar_a(x) (Tr[U_a(x) Ubar_a(x)] / N - 1)
        
        tr = one_ov_N * realtrace(&(s->link[mu]), &(s->link[mu])) - 1.0;
        tr *= dmu;
        scalar_mult_sum_adj_matrix(&(s->link[mu]), tr, &(s->f_U[mu])); // f_U[NUMLINK])
        
#endif
      }
    }
  }
  
 

  // Finally take adjoint and update the momentum
  // Include overall factor of kappa = N / (4lambda)
  // Subtract to reproduce -Adj(f_U)
  // Compute average gauge force in same loop
  
  //-edited------------
  
  tr = kappa * eps;  // eps = stepsize
  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      scalar_mult_dif_adj_matrix(&(s->f_U[mu]), tr, &(s->mom[mu]));
      returnit += realtrace(&(s->f_U[mu]), &(s->f_U[mu]));
    }
  }
  

  
  
  g_doublesum(&returnit);
  returnit *= kappa * kappa;

  // Add in force from separate determinant term if kappa_u1 non-zero
  if (kappa_u1 > IMAG_TOL)
    returnit += det_force(eps);

  return (eps * sqrt(returnit) / volume);
}
// -----------------------------------------------------------------



// Scalar force contribution phi



double scalar_force(Real eps) {



  register int i, mu, nu;
  register site *s;
  //char **local_pt[2][2];
  //---edited----- char **local_pt[2][2];
  char **local_pt[2][3]; 
  //int a, b, gather, flip = 0, index, next;
  double returnit = 0.0, tr;
  complex tcp;
  //matrix tmat, tmat2, *mat[2];
  
  matrix tmat,tmat2;
  
  msg_tag *tag[NUMLINK], *tag0[2], *tag1[2];
  
  
 

 // (1.0) contribution from 1/2*phi_commut^2, DmuUmu*phi_commut , [phi^bar,phi]*[varphi^bar,varphi],-2|[varphi,phi]|^2  wrt phi
 
                         // Zero force collectors
                         
  //----edited----                       
  FORALLSITES(i, s)
    clear_mat(&(s->f_phi));
     
   
    FORALLSITES(i, s) {
    
     
      mult_na(&(phi_commut[i]),  &(s->phi), &tmat);   // Initialize c= a*bdag => f_U_a[n] =phi_commut(n)phi^dag(n) // a = mu
      mult_an_dif(&(s->phi),&(phi_commut[i]),&tmat);  // c=c-adag*b => f_U_a[n] = phi_commut(n)phi^dag(n) - phi^dag(n)phi_commut(n)
      scalar_mult_sum_matrix(&tmat,1.0, &(s->f_phi)); // c = c + s * b => f_U[mu] = f_U[mu] + 
      
             
      mult_na(&(DmuUmu[i]),&(s->phi),&tmat);   // Initialize c= adag*b => f_U_a[n] = DmuUmu(n)phi^dag_a(n) // a = mu
      mult_an_dif(&(s->phi),&(DmuUmu[i]), &tmat); //c=c-a*bdag => f_U_a[n] = DmuUmu(n)phi^dag_a(n) - phi^dag_a(n)DmuUmu(n)
      scalar_mult_sum_matrix(&tmat,1.0, &(s->f_phi)); // c = c + s * b => f_U[mu] = f_U[mu]+
      
          
      mult_na(&(varphi_commut[i]),  &(s->phi), &tmat);   // Initialize c= a*bdag => f_U_a[n] =phi_commut(n)phi^dag(n) // a = mu
      mult_an_dif(&(s->phi),&(varphi_commut[i]),&tmat);  // c=c-adag*b => f_U_a[n] = phi_commut(n)phi^dag(n) - phi^dag(n)phi_commut(n)
      scalar_mult_sum_matrix(&tmat,1.0, &(s->f_phi)); // c = c + s * b => f_U[mu] = f_U[mu] + 
      
     
      mult_an(&(vp_commut[i]),&(s->varphi),&tmat);   // Initialize c= a*bdag => f_U_a[n] =phi_commut(n)phi^dag(n) // a = mu
      mult_na_dif(&(s->varphi),&(vp_commut[i]),&tmat);  // c=c-adag*b => f_U_a[n] = phi_commut(n)phi^dag(n) - phi^dag(n)phi_commut(n)
      scalar_mult_sum_matrix(&tmat,2.0, &(s->f_phi)); // c = c + s * b => f_U[mu] = f_U[mu]
      
     
     
    }
    

   
  
  
   
 // (2.0) contribution from 1/2*varphi_commut^2, DmuUmu*varphi_commut,  [phi^bar,phi]*[varphi^bar,varphi], -2|[varphi,phi]|^2 wrt varphi
 
                         // Zero force collectors
                         
  //----edited----                       
 FORALLSITES(i, s)
    clear_mat(&(s->f_varphi));
   

    FORALLSITES(i, s) {
    
     
        
      //clear_mat(&tmat);
      mult_na(&(varphi_commut[i]),  &(s->varphi), &tmat);   // Initialize c= a*bdag => f_U_a[n] =varphi_commut(n)varphi^dag(n) // a = mu
      mult_an_dif(&(s->varphi),&(varphi_commut[i]),&tmat);  // c=c-adag*b => f_U_a[n] = varphi_commut(n)varphi^dag(n) - varphi^dag(n)varphi_commut(n)
      scalar_mult_sum_matrix(&tmat,1.0, &(s->f_varphi)); // c = c + s * b => f_U[mu] = f_U[mu] + 
      
      
      mult_na(&(DmuUmu[i]),&(s->varphi),&tmat);   // Initialize c= adag*b => f_U_a[n] = DmuUmu(n)varphi^dag_a(n) // a = mu
      mult_an_dif(&(s->varphi),&(DmuUmu[i]), &tmat); //c=c-a*bdag => f_U_a[n] = DmuUmu(n)varphi^dag_a(n) - varphi^dag_a(n)DmuUmu(n)
      scalar_mult_sum_matrix(&tmat,1.0, &(s->f_varphi)); // c = c + s * b => f_U[mu] = f_U[mu]+
      
      
      mult_na(&(phi_commut[i]),  &(s->varphi), &tmat);   // Initialize c= a*bdag => f_U_a[n] =varphi_commut(n)varphi^dag(n) // a = mu
      mult_an_dif(&(s->varphi),&(phi_commut[i]),&tmat);  // c=c-adag*b => f_U_a[n] = varphi_commut(n)varphi^dag(n) - varphi^dag(n)varphi_commut(n)
      scalar_mult_sum_matrix(&tmat,1.0, &(s->f_varphi)); // c = c + s * b => f_U[mu] = f_U[mu] + 
      
      
     
      
      mult_na(&(s->phi),&(vp_commut[i]),&tmat);   // Initialize c= a*bdag => f_U_a[n] =phi_commut(n)phi^dag(n) // a = mu
      mult_an_dif(&(vp_commut[i]),&(s->phi),&tmat);  // c=c-adag*b => f_U_a[n] = phi_commut(n)phi^dag(n) - phi^dag(n)phi_commut(n)
      scalar_mult_sum_matrix(&tmat,2.0, &(s->f_varphi)); // c = c + s * b => f_U[mu] = f_U[mu] 
      
      
      
      } 
      
 
       
 
   // D+ original 

  //(3.0) contribution from -2|Dmuphi|^2 term wrt phi field
  
   for (mu = 0; mu < 3; mu++) 
   {
    local_pt[0][mu] = gen_pt[mu]; //  gen_pt[10 array element]
    local_pt[1][mu] = gen_pt[3 + mu];
   }

  FORALLDIR(mu) {
  
   tag[0] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                              goffset[mu]+1, EVENANDODD, local_pt[0][0] );// local_pt[0][0] = U_mu(n-mu)
   tag0[0] = start_gather_field(Dmuphi[mu], sizeof(matrix),
                              goffset[mu]+1, EVENANDODD, local_pt[1][0]);// local_pt[1][0] = DmuPhi[mu](n-mu)

    wait_gather(tag[0]);
    wait_gather(tag0[0]);
    
    FORALLSITES(i, s) {
      mult_an((matrix *)(local_pt[1][0][i]),(matrix *)(local_pt[0][0][i]), &tmat);   // Initialize c= adag*b => tmat = Dmuphi[mu]^dag(n-mu)U_mu(n-mu) // a = mu
      mult_na_dif(&(s->link[mu]),&(Dmuphi[mu][i]),  &tmat); //c=c-a*bdag => tmat = Dmuphi[mu]^dag(n-mu)U_mu(n-mu) - U_mu(n)Dmuphi[mu]^dag(n)
      scalar_mult_sum_matrix(&tmat,2.0, &(s->f_phi)); // c = c + s * b => f_phi = f_phi + Dmuphi[mu]^dag(n-mu)U_mu(n-mu) - U_mu(n)Dmuphi[mu]^dag(n)
    }
    
    
    cleanup_gather(tag[0]);
    cleanup_gather(tag0[0]);
  } 
  
  

   
   
       
  //(4.0) contribution from -2|Dmuvarphi|^2 term wrt varphi field
  
   for (mu = 0; mu < 3; mu++) 
   {
    local_pt[0][mu] = gen_pt[mu]; //  gen_pt[10 array element]
    local_pt[1][mu] = gen_pt[3 + mu];
   }

  FORALLDIR(mu) {
  
   tag[0] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                              goffset[mu]+1, EVENANDODD, local_pt[0][0] );// local_pt[0][0] = U_mu(n-mu)
   tag0[0] = start_gather_field(Dmuvarphi[mu], sizeof(matrix),
                              goffset[mu]+1, EVENANDODD, local_pt[1][0]);// local_pt[1][0] = DmuPhi[mu](n-mu)

    wait_gather(tag[0]);
    wait_gather(tag0[0]);
    
    FORALLSITES(i, s) {
      mult_an((matrix *)(local_pt[1][0][i]),(matrix *)(local_pt[0][0][i]), &tmat);   // Initialize c= adag*b => tmat = Dmuphi[mu]^dag(n-mu)U_mu(n-mu) // a = mu
      mult_na_dif(&(s->link[mu]),&(Dmuvarphi[mu][i]),  &tmat); //c=c-a*bdag => tmat = Dmuphi[mu]^dag(n-mu)U_mu(n-mu) - U_mu(n)Dmuphi[mu]^dag(n)
      scalar_mult_sum_matrix(&tmat,2.0, &(s->f_varphi)); // c = c + s * b => f_phi = f_phi + Dmuphi[mu]^dag(n-mu)U_mu(n-mu) - U_mu(n)Dmuphi[mu]^dag(n)
    }
    
    
    cleanup_gather(tag[0]);
    cleanup_gather(tag0[0]);
  } 
  
 
   

     // Only compute susy-breaking scalar potential term if bmass non-zero
  if (bmass > IMAG_TOL) {
    Real dmu,dmun;
#ifdef EIG_POT
    dmu = 2.0 * bmass * bmass;
    dmun = 1.0 * bmass * bmass;
    
#else
    Real trp,trv;
    dmu = 2.0 * one_ov_N * bmass * bmass;
   
   //dmu = 0.0; // edited
#endif

    FORALLSITES(i, s) {
     // FORALLDIR(mu) {
#ifdef EIG_POT
        
        
      
        // Ubar_a(x) [U_a(x) Ubar_a(x) - I]
        mult_na(&(s->phi), &(s->phi), &tmat);
        scalar_add_diag(&tmat, -1.0);
        scalar_mult_an_sum(&(s->phi), &tmat, dmu, &(s->f_phi));
        
        // Ubar_a(x) [U_a(x) Ubar_a(x) - I]
        mult_na(&(s->varphi), &(s->varphi), &tmat);
        scalar_add_diag(&tmat, -1.0);
        scalar_mult_an_sum(&(s->varphi), &tmat, dmu, &(s->f_varphi));
        
        

        
#else
        // Ubar_a(x) (Tr[U_a(x) Ubar_a(x)] / N - 1)
        trp =  one_ov_N * realtrace(&(s->phi), &(s->phi)) - 1.0;
        trp *= dmu;
        scalar_mult_sum_adj_matrix(&(s->phi), trp, &(s->f_phi)); // f_U[NUMLINK])
        trv =  one_ov_N * realtrace(&(s->varphi), &(s->varphi)) - 1.0;
        trv *= dmu;
        scalar_mult_sum_adj_matrix(&(s->varphi), trv, &(s->f_varphi)); // f_U[NUMLINK])
#endif
      //}
    }
  }
   
 

   
   
    
  // Finally take adjoint and update the momentum
  // Include overall factor of kappa = N / (4lambda)
  // Subtract to reproduce -Adj(f_U)
  // Compute average gauge force in same loop
  
  //-edited------------
  
  tr = kappa * eps;  // eps = stepsize
  

  
  
  FORALLSITES(i, s) {
    
      scalar_mult_dif_adj_matrix(&(s->f_phi), tr, &(s->mom_phi));                       //       c <-- c - s * bdag
      scalar_mult_dif_adj_matrix(&(s->f_varphi), tr, &(s->mom_varphi));               
         
          
                                     
      returnit += realtrace(&(s->f_phi), &(s->f_phi));                        //        
      returnit += realtrace(&(s->f_varphi), &(s->f_varphi)); 

    
  }
  
  g_doublesum(&returnit);
  returnit *= kappa * kappa;

  // Add in force from separate determinant term if kappa_u1 non-zero
  //if (kappa_u1 > IMAG_TOL)
    //returnit += det_force(eps);

  return (eps * sqrt(returnit) / volume);

}





// -----------------------------------------------------------------
// Separate routines for each term in the fermion force
// Plaquette determinant contributions to the fermion force
// Use Uinv, Udag_inv, UpsiU, Tr_Uinv and tr_dest for temporary storage
// Also use tempdet for temporary storage
// The accumulator names refer to the corresponding derivatives
// Assume compute_plaqdet() has already been run
// Appropriate adjoints set up in assemble_fermion_force
// A bit more code reuse may be possible
#ifdef SV
void detF(matrix *eta, matrix *psi[NUMLINK], int sign) {
  register int i;
  register site *s;
  int a, b, opp_b;
  Real localG = 0.5 * C2 * G;
  complex tc, tc2;
  msg_tag *mtag[8];
  matrix tmat;

  // Check sign while giving G proper sign
  if (sign == 1)
    localG *= -1.0;
  else if (sign != -1) {
    node0_printf("Error: incorrect sign in detF: %d\n", sign);
    terminate(1);
  }

  // Set up and store some basic ingredients
  FORALLSITES(i, s)
    tr_eta[i] = trace(&(eta[i]));
  // Need all directions for upcoming sums
  FORALLDIR(a) {
    // U_a(x)^{-1} and Udag_a(x)^{-1} are already in Uinv[a] and Udag_inv[a]
    // Save Tr[U_a(x)^{-1} psi_a(x) in Tr_Uinv[a]
    // Save U_a(x)^{-1} psi_a(x) U_a(x)^{-1} in UpsiU[a]
    FORALLSITES(i, s) {
      mult_nn(&(Uinv[a][i]), &(psi[a][i]), &tmat);
      mult_nn(&tmat, &(Uinv[a][i]), &(UpsiU[a][i]));
      Tr_Uinv[a][i] = trace(&tmat);

      // tempdet holds Tr[eta(x)] plaqdet[a][b](x)
      for (b = a + 1; b < NUMLINK; b++) {
        CMUL(tr_eta[i], plaqdet[a][b][i], tempdet[a][b][i]);
        CMUL(tr_eta[i], plaqdet[b][a][i], tempdet[b][a][i]);
      }
    }
  }

  // Now we are ready to gather, accumulate and add to force
  complex *plaq_term = malloc(sizeof *plaq_term * sites_on_node);
  complex *inv_term = malloc(sizeof *inv_term * sites_on_node);
  complex *adj_term = malloc(sizeof *adj_term * sites_on_node);

  // Now we are ready to gather, accumulate and add to force
  // TODO: Could try to overlap these gathers, but that looks nasty...
  FORALLDIR(a) {
    // Initialize accumulators for sums over b
    FORALLSITES(i, s) {
      plaq_term[i] = cmplx(0.0, 0.0);
      inv_term[i] = cmplx(0.0, 0.0);
      adj_term[i] = cmplx(0.0, 0.0);
    }
    FORALLDIR(b) {
      if (a == b)
        continue;

      // Summary of gathers and shorthand:
      //   D[a][b](x) is eta^D(x) tempdet[a][b](x)
      //   T[a](x) is Tr[U_a(x)^{-1} psi_a(x)]
      // 0) T[b](x + a - b) in two steps
      // 1) D[a][b](x - b)
      // 2) D[b][a](x - b)
      // 3) T[a](x + b)
      // 4) T[a](x - b)
      // 5) T[b](x + a)
      // 6) T[b](x - b)
      mtag[0] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[0]);
      mtag[1] = start_gather_field(tempdet[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[1]);
      mtag[2] = start_gather_field(tempdet[b][a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[2]);
      mtag[3] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                   goffset[b], EVENANDODD, gen_pt[3]);
      mtag[4] = start_gather_field(Tr_Uinv[a], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[4]);
      mtag[5] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[a], EVENANDODD, gen_pt[5]);
      mtag[6] = start_gather_field(Tr_Uinv[b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[6]);

      // Step two of Tr_Uinv[b](x - b + a) gather, including BC
      // Use tr_dest for temporary storage
      wait_gather(mtag[0]);
      FORALLSITES(i, s)
        CMULREAL(*((complex *)(gen_pt[0][i])), s->bc[a], tr_dest[i]);
      cleanup_gather(mtag[0]);
      mtag[0] = start_gather_field(tr_dest, sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

      // Now accumulate all three terms
      opp_b = OPP_LDIR(b);
      wait_gather(mtag[1]);         // 1) D[a][b](x - b)
      wait_gather(mtag[2]);         // 2) D[b][a](x - b)
      wait_gather(mtag[3]);         // 5) T[a](x + b)
      wait_gather(mtag[4]);         // 6) T[a](x - b)
      wait_gather(mtag[5]);         // 3) T[b](x + a)
      wait_gather(mtag[6]);         // 4) T[b](x - b)
      wait_gather(mtag[0]);         // 0) T[b](x + a - b)
      FORALLSITES(i, s) {
        // Accumulate plaq_term
        // D[b][a](x) {T[a](x) + T[b](x + a)}
        // gen_pt[5] is T[b](x + a)
        tc = *((complex *)(gen_pt[5][i]));
        tc2.real = Tr_Uinv[a][i].real + s->bc[a] * tc.real;
        tc2.imag = Tr_Uinv[a][i].imag + s->bc[a] * tc.imag;
        plaq_term[i].real += tempdet[b][a][i].real * tc2.real
                           - tempdet[b][a][i].imag * tc2.imag;
        plaq_term[i].imag += tempdet[b][a][i].imag * tc2.real
                           + tempdet[b][a][i].real * tc2.imag;

        // D[a][b](x - b) {T[a](x) + T[b](x - b)}
        // gen_pt[6] is T[b](x - b)
        tc = *((complex *)(gen_pt[6][i]));
        tc2.real = tc.real + s->bc[opp_b] * Tr_Uinv[a][i].real;
        tc2.imag = tc.imag + s->bc[opp_b] * Tr_Uinv[a][i].imag;
        // gen_pt[1] is D[a][b](x - b)
        tc = *((complex *)(gen_pt[1][i]));
        plaq_term[i].real += tc.real * tc2.real - tc.imag * tc2.imag;
        plaq_term[i].imag += tc.imag * tc2.real + tc.real * tc2.imag;

        // Accumulate adj_term
        // D[a][b](x) {T[b](x) + T[a](x + b)}
        // gen_pt[3] is T[a](x + b)
        tc = *((complex *)(gen_pt[3][i]));
        tc2.real = Tr_Uinv[b][i].real + s->bc[b] * tc.real;
        tc2.imag = Tr_Uinv[b][i].imag + s->bc[b] * tc.imag;
        adj_term[i].real += tempdet[a][b][i].real * tc2.real
                          - tempdet[a][b][i].imag * tc2.imag;
        adj_term[i].imag += tempdet[a][b][i].imag * tc2.real
                          + tempdet[a][b][i].real * tc2.imag;

        // D[b][a](x - b) {T[a](x - b) + T[b](x + a - b) bc[a](x - b)}
        // gen_pt[0] is T[b](x + a - b)
        // gen_pt[4] is T[a](x - b)
        CADD(*((complex *)(gen_pt[0][i])), *((complex *)(gen_pt[4][i])), tc);
        // gen_pt[2] is D[b][a](x - b)
        tc2 = *((complex *)(gen_pt[2][i]));
        adj_term[i].real += tc.real * tc2.real - tc.imag * tc2.imag;
        adj_term[i].imag += tc.imag * tc2.real + tc.real * tc2.imag;

        // Accumulate inv_term = sum_b D[b][a](x) + D[a][b](x - b)
        // gen_pt[1] is D[a][b](x - b)
        tc = *((complex *)(gen_pt[1][i]));
        inv_term[i].real += tempdet[b][a][i].real + s->bc[opp_b] * tc.real;
        inv_term[i].imag += tempdet[b][a][i].imag + s->bc[opp_b] * tc.imag;
      }
      cleanup_gather(mtag[0]);
      cleanup_gather(mtag[1]);
      cleanup_gather(mtag[2]);
      cleanup_gather(mtag[3]);
      cleanup_gather(mtag[4]);
      cleanup_gather(mtag[5]);
      cleanup_gather(mtag[6]);
    }

    // Now add to force
    FORALLSITES(i, s) {
      // Start with plaq_term hitting U_a(x)^{-1}
      CMULREAL(plaq_term[i], localG, tc);
      c_scalar_mult_sum_mat(&(Uinv[a][i]), &tc, &(s->f_U[a]));

      // Add adj_term hitting Udag_a(x)^{-1} followed by adjoint
      CMULREAL(adj_term[i], localG, tc);
      c_scalar_mult_sum_adj_mat(&(Udag_inv[a][i]), &tc, &(s->f_U[a]));

      // Finally subtract inv_term hitting U_a(x)^{-1} psi_a(x) U_a(x)^{-1}
      CMULREAL(inv_term[i], localG, tc);
      c_scalar_mult_dif_mat(&(UpsiU[a][i]), &tc, &(s->f_U[a]));
    }
  }
  free(plaq_term);
  free(inv_term);
  free(adj_term);
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Assemble fermion contributions to gauge link force,
//   f_U = Adj(Ms).D_U M(U, Ub).s - Adj[Adj(Ms).D_Ub M(U, Ub).s]
// "s" is sol while "Ms" is psol
// Copy these into persistent matrices for easier gathering
// Use tempmat, tempmat2, UpsiU, Tr_Uinv,
// tr_dest and Ddet[012] for temporary storage
// (many through calls to detF)
#ifndef PUREGAUGE
void assemble_fermion_force(Twist_Fermion *sol, Twist_Fermion *psol) {
  register int i, opp_mu, opp_rho;
  register site *s;
  char **local_pt[6][4]; // edited---- from local_pt[4][3]
  int mu, nu,rho, a, b, gather, flip = 0, index, next;

  msg_tag *mtag[NUMLINK], *tag0[6], *tag1[6], *tag2[6], *tag3[6];
 
  matrix *mat[2], tmat;
  
  for (int y = 0; y < 6; y++) {
    for (int x = 0; x < 4; x++) {
        local_pt[y][x] = gen_pt[y * 4 + x];
    }
  }
  


  mat[0] = tempmat;
  mat[1] = tempmat2;
  
   Real SGN;

  // For gathering it is convenient to copy the input Twist_Fermions
  // into persistent site, link and plaq fermions
  // We can reuse "src" and "dest" for this storage,
  // corresponding to "sol" and "psol", respectively
  FORALLSITES(i, s) {
    mat_copy(&(sol[i].Fsite), &(site_src[i]));
    adjoint(&(psol[i].Fsite), &(site_dest[i]));
    mat_copy(&(sol[i].Fsitez), &(sitez_src[i]));
    adjoint(&(psol[i].Fsitez), &(sitez_dest[i]));
    mat_copy(&(sol[i].Fsitezb), &(sitezb_src[i]));
    adjoint(&(psol[i].Fsitezb), &(sitezb_dest[i]));
    mat_copy(&(sol[i].Fsiteeb), &(siteeb_src[i]));
    adjoint(&(psol[i].Fsiteeb), &(siteeb_dest[i]));
    FORALLDIR(mu) {
      mat_copy(&(sol[i].Flink[mu]), &(link_src[mu][i]));
      adjoint(&(psol[i].Flink[mu]), &(link_dest[mu][i]));
      mat_copy(&(sol[i].Flinkb[mu]), &(linkb_src[mu][i]));
      adjoint(&(psol[i].Flinkb[mu]), &(linkb_dest[mu][i]));
      mat_copy(&(sol[i].Fthetab[mu]), &(thetab_src[mu][i]));
      adjoint(&(psol[i].Fthetab[mu]), &(thetab_dest[mu][i]));
    }
    for (mu = 0; mu < NPLAQ; mu++) {
      mat_copy(&(sol[i].Fplaq[mu]), &(plaq_src[mu][i]));
      adjoint(&(psol[i].Fplaq[mu]), &(plaq_dest[mu][i]));
    }
  }

#ifdef SV
  // Accumulate both terms in UpsiU[mu], use to initialize f_U[mu]
  // First calculate DUbar on eta Dbar_mu psi_mu (LtoS)
  // [psi_mu(x) eta(x + mu) - eta(x) psi_mu(x)]^dag
  mtag[0] = start_gather_field(site_dest, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_dest, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu], &tmat);
      mult_nn(&(link_src[mu][i]), &tmat, &(UpsiU[mu][i]));   // Initialize
      mult_nn_dif(&(site_dest[i]), &(link_src[mu][i]), &(UpsiU[mu][i]));
    }
    cleanup_gather(mtag[mu]);
  }

  // 2nd term, DUbar on psi_mu Dbar_mu eta (StoL)
  // [eta(x) psi_mu(x) - psi_mu(x) eta(x + mu)]^dag
  mtag[0] = start_gather_field(site_src, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]);
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(site_src, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu], &tmat);
      mult_nn_dif(&(link_dest[mu][i]), &tmat, &(UpsiU[mu][i]));
      mult_nn_sum(&(site_src[i]), &(link_dest[mu][i]), &(UpsiU[mu][i]));

      // Initialize the force collectors---done with UpsiU[mu]
 //edited    
      scalar_mult_adj_matrix(&(UpsiU[mu][i]), 0.5, &(s->f_U[mu]));
     // scalar_mult_sum_adj_matrix(&(UpsiU[mu][i]), -1.0, &(s->f_U[mu]));
    
    }
    cleanup_gather(mtag[mu]);
  }
  
  //(2)
  
  //term1
  //Sitezb to Lb
  
    mtag[0] = start_gather_field(sitezb_dest, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]);  //gen_pt[0]=sitezb_dest(n+0)
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(sitezb_dest, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);                   // gen_pt[mu]=sitezb_dest(n+mu)
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu], &tmat);    // tmat = s->bc[mu]*sitezb_dest(n+mu)
      mult_nn(&tmat,&(linkb_src[mu][i]), &(UpsiU[mu][i]));   // Initialize         UpsiU[mu][i] = s->bc[mu]zb(i+mu)*psi_mu(i)
      mult_nn_dif(&(linkb_src[mu][i]), &(sitezb_dest[i]), &(UpsiU[mu][i])); // UpsiU[mu][x] = s->bc[mu]zb(i+mu)*psi_mu(i) - psi_mu(i)zb(i)
    }
    cleanup_gather(mtag[mu]);
  }
  // Lb to Sitezb 
  // 2nd term,
  
  mtag[0] = start_gather_field(sitezb_src, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]); //gen_pt[0]=zb(n+0)
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(sitezb_src, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);                  //gen_pt[0]=zb(n+mu)
    }
    opp_mu = OPP_LDIR(mu);
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu]/*s->bc[opp_mu]*/, &tmat); // tmat = zb(n+mu)
      mult_nn_dif(&tmat, &(linkb_dest[mu][i]), &(UpsiU[mu][i])); //  UpsiU[mu][i] = UpsiU[mu][i] - zb(n+mu)psi_mu(n)
      mult_nn_sum(&(linkb_dest[mu][i]), &(sitezb_src[i]), &(UpsiU[mu][i])); // UpsiU[mu][i] = UpsiU[mu][i] - zb(n+mu)psi_mu(n) +psi*zb

      // Initialize the force collectors---done with UpsiU[mu]
     // scalar_mult_adj_matrix(&(UpsiU[mu][i]), 0.5, &(s->f_U[mu]));
    // scalar_mult_sum_adj_matrix(&(UpsiU[mu][i]), 0.5, &(s->f_U[mu]));         //b <-- b + s * a^dag
     scalar_mult_sum_matrix(&(UpsiU[mu][i]), 0.5, &(s->f_U[mu]));         //b <-- b + s * a
    
    }
    cleanup_gather(mtag[mu]); 
  }
  
  // (3)
  
  //Siteeb to Tb
  
    mtag[0] = start_gather_field(siteeb_dest, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]);  //gen_pt[0]=sitezb_dest(n+0)
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(siteeb_dest, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);                   // gen_pt[mu]=sitezb_dest(n+mu)
    }
    wait_gather(mtag[mu]);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu], &tmat);    // tmat = 1*sitezb_dest(n+mu)
      mult_nn(&tmat,&(thetab_src[mu][i]), &(UpsiU[mu][i]));   // Initialize         UpsiU[mu][i] =  zb(i+mu)*psi_mu(i)
      mult_nn_dif(&(thetab_src[mu][i]), &(siteeb_dest[i]), &(UpsiU[mu][i])); // UpsiU[mu][x] = zb(i+mu)*psi_mu(i) - psi_mu(i)zb(i)
    }
    cleanup_gather(mtag[mu]);
  }
  // Lb to Siteeb 
  // 2nd term,
  
  mtag[0] = start_gather_field(siteeb_src, sizeof(matrix),
                               goffset[0], EVENANDODD, gen_pt[0]); //gen_pt[0]=zb(n+0)
  FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {
      mtag[mu + 1] = start_gather_field(siteeb_src, sizeof(matrix),
                                        goffset[mu + 1], EVENANDODD,
                                        gen_pt[mu + 1]);                  //gen_pt[0]=zb(n+mu)
    }
    wait_gather(mtag[mu]);
    opp_mu = OPP_LDIR(mu);
    FORALLSITES(i, s) {
      scalar_mult_matrix((matrix *)(gen_pt[mu][i]), s->bc[mu]/*s->bc[opp_mu]*/, &tmat); // tmat = zb(n+mu)
      mult_nn_dif(&tmat, &(thetab_dest[mu][i]), &(UpsiU[mu][i])); //  UpsiU[mu][i] = UpsiU[mu][i] - zb(n+mu)psi_mu(n)
      mult_nn_sum(&(thetab_dest[mu][i]), &(siteeb_src[i]), &(UpsiU[mu][i])); // UpsiU[mu][i] = UpsiU[mu][i] - zb(n+mu)psi_mu(n) +psi*zb

      scalar_mult_sum_matrix(&(UpsiU[mu][i]), 0.5, &(s->f_U[mu]));         //b <-- b + s * a
    
    }
    cleanup_gather(mtag[mu]); 
  }
  
 
  
#else
  FORALLDIR(mu) {                         // Zero force collectors
    FORALLSITES(i, s)
      clear_mat(&(s->f_U[mu]));
  }
#endif
#ifdef VP
  // Now calculate DU on chi_{munu} D_mu(U) psi_nu
  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(link_src[1], sizeof(matrix),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  // Prepare and gather other term in tempmat*
  index = plaq_index[0][1];
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_nn(&(plaq_dest[index][i]), &(link_src[1][i]), &(mat[0][i]));
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) {   // Start next gathers
        if (nu == NUMLINK - 1) {
          a = mu + 1;
          b = 0;
        }
        else if (nu == mu - 1) {
          a = mu;
          b = nu + 2;
        }
        else {
          a = mu;
          b = nu + 1;
        }
        tag0[gather] = start_gather_field(link_src[b], sizeof(matrix),
                                          goffset[a], EVENANDODD,
                                          local_pt[gather][0]);

        // Prepare and gather other term in tempmat*
        next = plaq_index[a][b];
        FORALLSITES(i, s) {
          if (a > b) {      // plaq_dest is anti-symmetric under a <--> b
            scalar_mult_matrix(&(plaq_dest[next][i]), -1.0, &tmat);
          }                 // Suppress compiler error
          else
            mat_copy(&(plaq_dest[next][i]), &tmat);
          mult_nn(&tmat, &(link_src[b][i]), &(mat[gather][i]));
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {    // plaq_dest is anti-symmetric under mu <--> nu
          scalar_mult_matrix((matrix *)(local_pt[flip][0][i]),
                             s->bc[mu], &tmat);
        }                 // Suppress compiler error
        else
          scalar_mult_matrix((matrix *)(local_pt[flip][0][i]),
                             -1.0 * s->bc[mu], &tmat);

        mult_nn_sum(&tmat, &(plaq_dest[index][i]), &(s->f_U[mu]));
        sum_matrix((matrix *)(local_pt[flip][1][i]), &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }

  // 2nd term
  // Start first set of gathers (mu = 0 and nu = 1)
  tag0[0] = start_gather_field(link_dest[1], sizeof(matrix),
                               goffset[0], EVENANDODD, local_pt[0][0]);

  // Prepare and gather other term in tempmat*
  index = plaq_index[0][1];
  FORALLSITES(i, s)     // mu = 0 < nu = 1
    mult_nn(&(plaq_src[index][i]), &(link_dest[1][i]), &(mat[0][i]));
  tag1[0] = start_gather_field(mat[0], sizeof(matrix),
                               goffset[1] + 1, EVENANDODD, local_pt[0][1]);

  // Main loop
  FORALLDIR(mu) {
    FORALLDIR(nu) {
      if (mu == nu)
        continue;

      gather = (flip + 1) % 2;
      if (mu < NUMLINK - 1 || nu < NUMLINK - 2) {   // Start next gathers
        if (nu == NUMLINK - 1) {
          a = mu + 1;
          b = 0;
        }
        else if (nu == mu - 1) {
          a = mu;
          b = nu + 2;
        }
        else {
          a = mu;
          b = nu + 1;
        }
        tag0[gather] = start_gather_field(link_dest[b], sizeof(matrix),
                                          goffset[a], EVENANDODD,
                                          local_pt[gather][0]);

        // Prepare and gather other term in tempmat*
        next = plaq_index[a][b];
        FORALLSITES(i, s) {
          if (a > b) {      // plaq_src is anti-symmetric under a <--> b
            scalar_mult_matrix(&(plaq_src[next][i]), -1.0, &tmat);
          }                 // Suppress compiler error
          else
            mat_copy(&(plaq_src[next][i]), &tmat);
          mult_nn(&tmat, &(link_dest[b][i]), &(mat[gather][i]));
        }
        tag1[gather] = start_gather_field(mat[gather], sizeof(matrix),
                                          goffset[b] + 1, EVENANDODD,
                                          local_pt[gather][1]);
      }

      index = plaq_index[mu][nu];
      wait_gather(tag0[flip]);
      wait_gather(tag1[flip]);
      FORALLSITES(i, s) {
        if (mu > nu) {    // plaq_src is anti-symmetric under mu <--> nu
          scalar_mult_matrix(&(plaq_src[index][i]), -1.0 * s->bc[mu], &tmat);
        }                 // Suppress compiler error
        else
          scalar_mult_matrix(&(plaq_src[index][i]), s->bc[mu], &tmat);

        mult_nn_sum((matrix *)(local_pt[flip][0][i]), &tmat, &(s->f_U[mu]));
        dif_matrix((matrix *)(local_pt[flip][1][i]), &(s->f_U[mu]));
      }
      cleanup_gather(tag0[flip]);
      cleanup_gather(tag1[flip]);
      flip = gather;
    }
  }
  
  
    //(2) chi->linkb
  
    // Term1 PtoLb
    
  Real SIGN_array[3] = {-1.0, 1.0, -1.0};  
    
  FORALLSITES(i, s)
    clear_mat(&(s->f_phi));  
   

  
  FORALLDIR(mu) {
  
  
       Real SIGN = SIGN_array[mu];
        

        tag0[mu] = start_gather_field(plaq_src[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][mu]);           // chi_mu(n+mu)

 
        tag1[mu] = start_gather_field(linkb_dest[mu], sizeof(matrix),
                                          goffset[mu]+1, EVENANDODD,
                                          local_pt[1][mu]);           // Lb_mu(n-mu)
      

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      opp_mu = OPP_LDIR(mu);
      
      FORALLSITES(i, s) {
        
       mult_nn((matrix *)(local_pt[0][mu][i]),&(linkb_dest[mu][i]), &(UpsiU[mu][i]) );          // &(UpsiU[mu][i])  = chi_mu(n+mu)Lb(n)
       scalar_mult_matrix(&(UpsiU[mu][i]), s->bc[mu], &(UpsiU[mu][i]));                         // &(UpsiU[mu][i])  = bc[mu]chi_mu(n+mu)Lb(n)
       scalar_mult_nn_dif((matrix *)(local_pt[1][mu][i]),&(plaq_src[mu][i]),s->bc[opp_mu],&(UpsiU[mu][i]) );      // &(UpsiU[mu][i])  = bc[mu]chi_mu(n+mu)Lb(n) - bc[opp_mu]Lb_mu(n-mu)chi_mu(n)

       
         
      }
      
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
     
    
  }

  // 2nd term LbtoP


  FORALLDIR(mu) {
  
  
       Real SIGN = SIGN_array[mu];
        

        tag0[mu] = start_gather_field(plaq_dest[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][mu]);           // chi_mu(n+mu)

 
        tag1[mu] = start_gather_field(linkb_src[mu], sizeof(matrix),
                                          goffset[mu]+1, EVENANDODD,
                                          local_pt[1][mu]);           // Lb(n-mu)
      

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      opp_mu = OPP_LDIR(mu);
      
      FORALLSITES(i, s) {
        
       scalar_mult_nn_sum((matrix *)(local_pt[1][mu][i]),&(plaq_dest[mu][i]),s->bc[opp_mu],&(UpsiU[mu][i]) );          // &(UpsiU[mu][i])  += bc[opp_mu]Lb(n-mu)chi_mu(n)
       scalar_mult_nn_dif((matrix *)(local_pt[0][mu][i]),&(linkb_src[mu][i]),s->bc[mu],&(UpsiU[mu][i]) );      // &(UpsiU[mu][i])  = bc[opp_mu]Lb(n-mu)chi_mu(n) - bc[mu]chi_mu(n+mu)Lb_mu(n)
       scalar_mult_matrix(&(UpsiU[mu][i]), SIGN, &(UpsiU[mu][i]) );                       // &(UpsiU[mu][i]) += SIGN[mu} * (bc[opp_mu]Lb(n-mu)chi_mu(n) - bc[mu]chi_mu(n+mu)Lb_mu(n))

       
       scalar_mult_sum_adj_matrix(&(UpsiU[mu][i]), -1.0, &(s->f_phi));
       
         
      }
      
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
     
    
  }
  
  
      //(3) chi->Thetab
  
    // Term1 PtoTb
    
  FORALLSITES(i, s)
    clear_mat(&(s->f_varphi));  
   

  
  FORALLDIR(mu) {
  
  
       Real SIGN = SIGN_array[mu];
        

        tag0[mu] = start_gather_field(plaq_src[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][mu]);           // chi_mu(n+mu)

 
        tag1[mu] = start_gather_field(thetab_dest[mu], sizeof(matrix),
                                          goffset[mu]+1, EVENANDODD,
                                          local_pt[1][mu]);           // Lb_mu(n-mu)
      

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      opp_mu = OPP_LDIR(mu);
      
      FORALLSITES(i, s) {
        
       mult_nn((matrix *)(local_pt[0][mu][i]),&(thetab_dest[mu][i]), &(UpsiU[mu][i]) );      // &(UpsiU[mu][i])  = chi_mu(n+mu)Lb(n)
       scalar_mult_matrix(&(UpsiU[mu][i]), s->bc[mu], &(UpsiU[mu][i]));                      // &(UpsiU[mu][i])  = bc[mu]chi_mu(n+mu)Lb(n)
       scalar_mult_nn_dif((matrix *)(local_pt[1][mu][i]),&(plaq_src[mu][i]),s->bc[opp_mu],&(UpsiU[mu][i]) );      // &(UpsiU[mu][i])  = bc[mu]chi_mu(n+mu)Lb(n) - bc[op_mu]Lb_mu(n-mu)chi_mu(n)

       
         
      }
      
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
     
    
  }

  // 2nd term TbtoP


  FORALLDIR(mu) {
  
  
       Real SIGN = SIGN_array[mu];
        

        tag0[mu] = start_gather_field(plaq_dest[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][mu]);           // chi_mu(n+mu)

 
        tag1[mu] = start_gather_field(thetab_src[mu], sizeof(matrix),
                                          goffset[mu]+1, EVENANDODD,
                                          local_pt[1][mu]);           // Lb(n-mu)
      

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      opp_mu = OPP_LDIR(mu);
      
      FORALLSITES(i, s) {
        
       scalar_mult_nn_sum((matrix *)(local_pt[1][mu][i]),&(plaq_dest[mu][i]),s->bc[opp_mu],&(UpsiU[mu][i]) );          // &(UpsiU[mu][i])  += bc[opp_mu]Lb(n-mu)chi_mu(n)
       scalar_mult_nn_dif((matrix *)(local_pt[0][mu][i]),&(thetab_src[mu][i]),s->bc[mu],&(UpsiU[mu][i]) );      // &(UpsiU[mu][i])  = bc[opp_mu]Lb(n-mu)chi_mu(n) - bc[mu]chi_mu(n+mu)Lb_mu(n)
       scalar_mult_matrix(&(UpsiU[mu][i]), -1.0*SIGN, &(UpsiU[mu][i]) );                       //&(UpsiU[mu][i]) += SIGN[mu} * (bc[opp_mu]Lb(n-mu)chi_mu(n) - bc[mu]chi_mu(n+mu)Lb_mu(n))

       
       scalar_mult_sum_adj_matrix(&(UpsiU[mu][i]), -1.0, &(s->f_varphi));
       
         
      }
      
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
     
    
  }
  
  

  
#endif

#ifdef SV
  // Plaquette determinant contributions if G is non-zero
  if (doG) {
    // First connect link_src with site_dest[DIMF - 1]^dag (LtoS)
    detF(site_dest, link_src, PLUS);

    // Second connect site_src[DIMF - 1] with link_dest^dag (StoL)
    detF(site_src, link_dest, MINUS);
  }
#endif



#ifdef SS


  // T1

//(1)  
  // First calculate SztoSeb 
  // [zeta -> etab  ]

  FORALLSITES(i, s) {
   
     
   
      mult_nn(&(sitez_src[i]),&(siteeb_dest[i]),&(tempmat[i]));   
      mult_nn_dif(&(siteeb_dest[i]),&(sitez_src[i]),&(tempmat[i]));  
      
     // scalar_mult_matrix(&tmat, -1.0, &tmat); 
      
     // 2nd calculate SebtoSz 
     // [etab -> zeta ] 
     
      
      mult_nn_sum(&(siteeb_src[i]), &(sitez_dest[i]), &(tempmat[i]));   
      mult_nn_dif(&(sitez_dest[i]), &(siteeb_src[i]), &(tempmat[i]));     
      scalar_mult_sum_matrix(&(tempmat[i]), -1.0, &(s->f_varphi));  // s->f_varphi = f_varphi(i) at s/i site structure 
      }
  
  //(4)  
  // First calculate SztoSeb 
  // [eta -> zetab  ]

   FORALLSITES(i, s) {
   
      mult_nn(&(site_src[i]),&(sitezb_dest[i]),&tmat);   
      mult_nn_dif(&(sitezb_dest[i]),&(site_src[i]),&tmat);  
      
     // scalar_mult_matrix(&tmat, -1.0, &tmat); 
      
     // 2nd calculate SebtoSz 
     // [zetab -> eta ] 
     
      
      mult_nn_sum(&(sitezb_src[i]), &(site_dest[i]), &tmat);   
      mult_nn_dif(&(site_dest[i]), &(sitezb_src[i]), &tmat);     
      scalar_mult_sum_adj_matrix(&tmat, -0.5, &(s->f_varphi));  // s->f_varphi = f_varphi(i) at s/i site structure 
      
      
      
    }   
    
    
    // T2


  
//(2)  
  // First calculate SztoSeb 
  // [zeta -> zetab  ]

   FORALLSITES(i, s) {
   
      mult_nn(&(sitez_src[i]),&(sitezb_dest[i]),&tmat);   
      mult_nn_dif(&(sitezb_dest[i]),&(sitez_src[i]),&tmat);  
      

      
     // 2nd calculate SebtoSz 
     // [zetab -> zeta ] 
     
      
      mult_nn_sum(&(sitezb_src[i]), &(sitez_dest[i]), &tmat);   
      mult_nn_dif(&(sitez_dest[i]), &(sitezb_src[i]), &tmat);     
      scalar_mult_sum_matrix(&tmat, 1.0, &(s->f_phi));  // s->f_varphi = f_varphi(i) at s/i site structure 
      
      
      
    }
 
//(3)  
  // First calculate SztoSeb 
  // [eta -> etab  ]

   FORALLSITES(i, s) {
   
      mult_nn(&(site_src[i]),&(siteeb_dest[i]),&tmat);   
      mult_nn_dif(&(siteeb_dest[i]),&(site_src[i]),&tmat);  
      
 
      
     // 2nd calculate SebtoSz 
     // [etab -> eta ] 
     
      
      mult_nn_sum(&(siteeb_src[i]), &(site_dest[i]), &tmat);   
      mult_nn_dif(&(site_dest[i]), &(siteeb_src[i]), &tmat);     
      scalar_mult_sum_adj_matrix(&tmat, -0.5, &(s->f_phi));  // s->f_varphi = f_varphi(i) at s/i site structure 
      
      
      
    } 
       
   
   
 
 
#else
                       // Zero force collectors
    FORALLSITES(i, s)
      clear_mat(&(s->f_varphi));
    FORALLSITES(i, s)
      clear_mat(&(s->f_phi));  
 
      
#endif


#ifdef VV

  
 
  // T1

  
//(1)  
  // 1st calculate LtoLb 
  // [psi -> psib  ]
  

  FORALLSITES(i, s) {           // Set up first gather
    mult_nn(&(linkb_dest[0][i]),&(link_src[0][i]), &(mat[0][i])); //mat[mu][i]=psib[mu]*psi[mu]
    }
  
  mtag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[mu]=psib[mu](n-mu)*psi[mu](n-mu)
  
    FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_nn(&(linkb_dest[nu][i]),&(link_src[nu][i]), &(mat[gather][i])); //mat[mu][i]=psib[mu]*psi[mu]
      mtag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); //gen_pt[mu]=psib[mu](n-mu)*psi[mu](n-mu)
    }

    wait_gather(mtag[mu]);
    opp_mu = OPP_LDIR(mu);
    
    FORALLSITES(i, s) {
    
      mult_nn(&(link_src[mu][i]),&(linkb_dest[mu][i]), &(UpsiU[mu][i]));                                // UpsiU[mu] = psi_mu(n)psib_mu(n)
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]),s->bc[opp_mu]*s->bc[opp_mu], &(UpsiU[mu][i])); // UpsiU[mu] = psi_mu(n)psib_mu(n) - bc[opp_mu]*bc[opp_mu]*psib[mu](n-mu)*psi[mu](n-mu)
    }
    
    cleanup_gather(mtag[mu]);
    flip = gather;
  }  
  
  // 2nd calculate LbtoL 
  // [psib -> psi  ]

   

  FORALLSITES(i, s) {           // Set up first gather
    mult_nn(&(linkb_src[0][i]),&(link_dest[0][i]), &(mat[0][i])); //mat[mu][i]=psib[mu]*psi[mu]
    }
  
  mtag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[mu]=psib[mu](n-mu)*psi[mu](n-mu)
  
    FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_nn(&(linkb_src[nu][i]),&(link_dest[nu][i]), &(mat[gather][i])); //mat[mu][i]=psib[mu]*psi[mu]
      mtag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); //gen_pt[mu]=psib[mu](n-mu)*psi[mu](n-mu)
    }

    wait_gather(mtag[mu]);
    opp_mu = OPP_LDIR(mu);
    FORALLSITES(i, s) {
    
      mult_nn_dif(&(link_dest[mu][i]),&(linkb_src[mu][i]), &(UpsiU[mu][i]));     // UpsiU[mu] = UpsiU - psi_mu(n)psib_mu(n)
      scalar_mult_sum_matrix((matrix *)(gen_pt[mu][i]),s->bc[opp_mu]*s->bc[opp_mu], &(UpsiU[mu][i]));// UpsiU[mu] = UpsiU-psi_mu(n)psib_mu(n)+bc[opp_mu]*bc[opp_mu]*psib[mu](n-mu)*psi[mu](n-mu)
     
      scalar_mult_sum_matrix(&(UpsiU[mu][i]), 1.0, &(s->f_varphi));         //b <-- b + s * a
    }
    
    cleanup_gather(mtag[mu]);
    flip = gather;
  }
  
 
  
  
  
  //--------------To be edited by replacing Lb with Tb-----------------------------
  
  //(2)  
  
  // 1st calculate LtoTb 
  // [psi -> thetab  ]
  

  FORALLSITES(i, s) {           // Set up first gather
    mult_nn(&(thetab_dest[0][i]),&(link_src[0][i]), &(mat[0][i])); //mat[mu][i]=thetab[mu]*psi[mu]
    }
  
  mtag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[mu]=thetab[mu](n-mu)*psi[mu](n-mu)
  
    FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_nn(&(thetab_dest[nu][i]),&(link_src[nu][i]), &(mat[gather][i])); //mat[mu][i]=thetab[mu]*psi[mu]
      mtag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); //gen_pt[mu]=thetab[mu](n-mu)*psi[mu](n-mu)
    }

    wait_gather(mtag[mu]);
    opp_mu = OPP_LDIR(mu);
    FORALLSITES(i, s) {
    
      mult_nn(&(link_src[mu][i]),&(thetab_dest[mu][i]), &(UpsiU[mu][i]));   // UpsiU[mu] = psi_mu(n)thetab_mu(n)
      scalar_mult_dif_matrix((matrix *)(gen_pt[mu][i]),s->bc[opp_mu]*s->bc[opp_mu], &(UpsiU[mu][i])); // UpsiU[mu] = psi_mu(n)thetab_mu(n) - bc[opp_mu]*bc[opp_mu]*thetab[mu](n-mu)*psi[mu](n-mu)
    }
    
    cleanup_gather(mtag[mu]);
    flip = gather;
  }  
  
  // 2nd calculate TbtoL 
  
  // [Thetab -> psi  ]

   

  FORALLSITES(i, s) {           // Set up first gather
    mult_nn(&(thetab_src[0][i]),&(link_dest[0][i]), &(mat[0][i])); //mat[mu][i]=thetab[mu]*psi[mu]
    }
  
  mtag[0] = start_gather_field(mat[0], sizeof(matrix),
                              goffset[0] + 1, EVENANDODD, gen_pt[0]);  //gen_pt[mu]=thetab[mu](n-mu)*psi[mu](n-mu)
  
    FORALLDIR(mu) {
    if (mu < NUMLINK - 1) {   // Start next gather
      nu = mu + 1;
      gather = (flip + 1) % 2;
      FORALLSITES(i, s)
        mult_nn(&(thetab_src[nu][i]),&(link_dest[nu][i]), &(mat[gather][i])); //mat[mu][i]=thetab[mu]*psi[mu]
      mtag[nu] = start_gather_field(mat[gather], sizeof(matrix),
                                   goffset[nu] + 1, EVENANDODD, gen_pt[nu]); //gen_pt[mu]=thetab[mu](n-mu)*psi[mu](n-mu)
    }

    wait_gather(mtag[mu]);
    opp_mu = OPP_LDIR(mu);
    FORALLSITES(i, s) {
    
      mult_nn_dif(&(link_dest[mu][i]),&(thetab_src[mu][i]), &(UpsiU[mu][i]));     // UpsiU[mu] = UpsiU - psi_mu(n)thetab_mu(n)
      scalar_mult_sum_matrix((matrix *)(gen_pt[mu][i]),s->bc[opp_mu]*s->bc[opp_mu], &(UpsiU[mu][i]));// UpsiU[mu]=UpsiU-psi_mu(n)thetab_mu(n) + bc[opp_mu]*bc[opp_mu]thetab[mu](n-mu)*psi[mu](n-mu)
     
      scalar_mult_sum_matrix(&(UpsiU[mu][i]), 1.0, &(s->f_phi));         //b <-- b + s * a
    }
    
    cleanup_gather(mtag[mu]);
    flip = gather;
  }
  

  
    //(3) psib-D-thetab
  
  
// 1st calculate thetab->linkb     
  
//--------  
  FORALLDIR(mu) {                         // Zero force collectors
    FORALLSITES(i, s)
      clear_mat(&(UpsiU[mu][i]));
  }
 

  int mu_vals[6]  = {0, 0, 1, 1, 2, 2};
  int nu_vals[6]  = {1, 2, 2, 0, 0, 1};
  int rho_vals[6] = {2, 1, 0, 2, 1, 0}; 
  Real sgn_val_dp[6] = {-1,1,-1,1,-1,1};
  

     mu  = mu_vals[0];
     nu  = nu_vals[0];
     rho = rho_vals[0];
     
          

        tag0[0] = start_gather_field(thetab_src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][0]);              // thetab_rho(n-rho)      
                                          
        tag1[0] = start_gather_field(thetab_src[rho], sizeof(matrix),
                                          goffset[nu], EVENANDODD,
                                          local_pt[0][1]);              // thetab_rho(n+nu)      
                                  

        tag2[0] = start_gather_field(linkb_dest[mu], sizeof(matrix),  
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[0][2]);                // psib_mu(n-mu)
                                         
        tag3[0] = start_gather_field(linkb_dest[mu], sizeof(matrix),
                                         goffset[nu], EVENANDODD,
                                         local_pt[0][3]);                // psib_mu(n+nu)   
 

  for (int j = 0; j < 6; j++) {
  
  
  if(j<5) 
    
  {
     mu  = mu_vals[j+1];
     nu  = nu_vals[j+1];
     rho = rho_vals[j+1];
     
          

        tag0[j+1] = start_gather_field(thetab_src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[j+1][0]);              // thetab_rho(n-rho)      
                                          
        tag1[j+1] = start_gather_field(thetab_src[rho], sizeof(matrix),
                                          goffset[nu], EVENANDODD,
                                          local_pt[j+1][1]);              // thetab_rho(n+nu)      
                                  

        tag2[j+1] = start_gather_field(linkb_dest[mu], sizeof(matrix),  
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[j+1][2]);                // psib_mu(n-mu)
                                         
        tag3[j+1] = start_gather_field(linkb_dest[mu], sizeof(matrix),
                                         goffset[nu], EVENANDODD,
                                         local_pt[j+1][3]);                // psib_mu(n+nu)         
                                         
   }
   
    
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];   
     SGN = sgn_val_dp[j];                                      
                                                                          
      wait_gather(tag0[j]);
      wait_gather(tag1[j]);
      wait_gather(tag2[j]);
      wait_gather(tag3[j]);
      
      opp_rho = OPP_LDIR(rho);
      opp_mu = OPP_LDIR(mu);
      
      FORALLSITES(i, s) {
       
       mult_nn((matrix *)(local_pt[j][2][i]),(matrix *)(local_pt[j][1][i]), &tmat);             // UpsiU = psib_mu(n-mu)*thetab_rho(n+nu)
       scalar_mult_matrix(&tmat, s->bc[opp_mu], &tmat);                                         // UpsiU = bc[opp_mu]psib_mu(n-mu)*thetab_rho(n+nu)
       scalar_mult_matrix(&tmat, s->bc[nu], &tmat);                                             // UpsiU = bc[opp_mu]bc[nu]*psib_mu(n-mu)*thetab_rho(n+nu)
       scalar_mult_nn_dif((matrix *)(local_pt[j][0][i]),(matrix *)(local_pt[j][3][i]),s->bc[opp_rho] * s->bc[nu], &tmat);
      
       
       scalar_mult_sum_matrix(&tmat,1.0*SGN, &(UpsiU[nu][i])); // dest[i] = dest[i] - 1.0*{psib_mu(n-mu)*thetab_rho(n+nu) -thetab_rho(n-rho)psib_mu(n+nu)}
    

      }
      
      cleanup_gather(tag0[j]);
      cleanup_gather(tag1[j]);
      cleanup_gather(tag2[j]);
      cleanup_gather(tag3[j]);
  
    
  
      
    }
    
// thetab - D - psib    

// 2nd calculate linkb -> thetab     
  

  
  Real sgn_val_dm[6] = {1,-1, 1,-1, 1,-1}; 
  

     mu  = mu_vals[0];
     nu  = nu_vals[0];
     rho = rho_vals[0];
     
          

        tag0[0] = start_gather_field(linkb_src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[0][0]);              // linkb_rho(n-rho)      
                                          
        tag1[0] = start_gather_field(linkb_src[rho], sizeof(matrix),
                                          goffset[nu], EVENANDODD,
                                          local_pt[0][1]);              // linkb_rho(n+nu)      
                                  

        tag2[0] = start_gather_field(thetab_dest[mu], sizeof(matrix),  
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[0][2]);                // thetab_mu(n-mu)
                                         
        tag3[0] = start_gather_field(thetab_dest[mu], sizeof(matrix),
                                         goffset[nu], EVENANDODD,
                                         local_pt[0][3]);                // thetab_mu(n+nu)         
 

  for (int j = 0; j < 6; j++) {
 
   if(j<5)  
   
  {
     mu  = mu_vals[j+1];
     nu  = nu_vals[j+1];
     rho = rho_vals[j+1];
     
          

        tag0[j+1] = start_gather_field(linkb_src[rho], sizeof(matrix),
                                          goffset[rho]+1, EVENANDODD,
                                          local_pt[j+1][0]);              // linkb_rho(n-rho)      
                                          
        tag1[j+1] = start_gather_field(linkb_src[rho], sizeof(matrix),
                                          goffset[nu], EVENANDODD,
                                          local_pt[j+1][1]);              // linkb_rho(n+nu)      
                                  

        tag2[j+1] = start_gather_field(thetab_dest[mu], sizeof(matrix),  
                                         goffset[mu] + 1, EVENANDODD,
                                         local_pt[j+1][2]);                // thetab_mu(n-mu)
                                         
        tag3[j+1] = start_gather_field(thetab_dest[mu], sizeof(matrix),
                                         goffset[nu], EVENANDODD,
                                         local_pt[j+1][3]);                // thetab_mu(n+nu)     
                                         
   }

   
                                     
     mu  = mu_vals[j];
     nu  = nu_vals[j];
     rho = rho_vals[j];   
     SGN = sgn_val_dm[j];                                              
                                                                          
      wait_gather(tag0[j]);
      wait_gather(tag1[j]);
      wait_gather(tag2[j]);
      wait_gather(tag3[j]);
      
      opp_rho = OPP_LDIR(rho);
      opp_mu = OPP_LDIR(mu);
      
      FORALLSITES(i, s) {
       
       mult_nn((matrix *)(local_pt[j][2][i]),(matrix *)(local_pt[j][1][i]), &tmat);          // UpsiU = thetab_mu(n-mu)*Lb_rho(n+nu)
       scalar_mult_matrix(&tmat, s->bc[opp_mu], &tmat);                                      // UpsiU = bc[opp_mu]thetab_mu(n-mu)*Lb_rho(n+nu)
       scalar_mult_matrix(&tmat, s->bc[nu], &tmat);                                          // UpsiU = bc[opp_mu]bc[nu]thetab_mu(n-mu)*Lb_rho(n+nu)
       scalar_mult_nn_dif((matrix *)(local_pt[j][0][i]),(matrix *)(local_pt[j][3][i]),s->bc[opp_rho] * s->bc[nu],&tmat); 
       
       
       
       scalar_mult_sum_matrix(&tmat,1.0*SGN, &(UpsiU[nu][i])); // dest[i] = dest[i] - 1.0*{thetab_mu(n-mu)*Lb_rho(n+nu) -Lb_rho(n-rho)thetab_mu(n+nu)}
      

      }
      
     
      
      cleanup_gather(tag0[j]);
      cleanup_gather(tag1[j]);
      cleanup_gather(tag2[j]);
      cleanup_gather(tag3[j]);
  
     
  
      
    }  
  
  
    FORALLDIR(mu) { 
       FORALLSITES(i, s) {
        scalar_mult_sum_adj_matrix(&(UpsiU[mu][i]), -1.0, &(s->f_U[mu]));
                       }
    }

  
    
 #else

#endif


#ifdef SP
 


  // Term1 PtoSz
  
  FORALLDIR(mu) {
  
  
     //  Real SIGN = SIGN_array[mu];
        

        tag0[mu] = start_gather_field(plaq_src[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][mu]);           // chi_mu(n+mu)

 
        tag1[mu] = start_gather_field(sitez_dest, sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);           // Z(n+mu)
      

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      
      FORALLSITES(i, s) {
        
       mult_nn(&(sitez_dest[i]),(matrix *)(local_pt[0][mu][i]), &(UpsiU[mu][i]) );          // &(UpsiU[mu][i])  = Z(n)chi_mu(n+mu)
       scalar_mult_matrix(&(UpsiU[mu][i]) , s->bc[mu], &(UpsiU[mu][i]) );                   // &(UpsiU[mu][i])  = bc[mu]Z(n)chi_mu(n+mu)
       scalar_mult_nn_dif((matrix *)(local_pt[0][mu][i]), (matrix *)(local_pt[1][mu][i]),s->bc[mu] * s->bc[mu], &(UpsiU[mu][i]) ); 
 
         
      }
      
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
     
    
  }

  // 2nd term StoP


  FORALLDIR(mu) {
  
  
        Real SIGN = SIGN_array[mu];
        

        tag0[mu] = start_gather_field(plaq_dest[mu], sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[0][mu]);           // chi_mu(n+mu)

 
        tag1[mu] = start_gather_field(sitez_src, sizeof(matrix),
                                          goffset[mu], EVENANDODD,
                                          local_pt[1][mu]);           // Z(n+mu)
      

      wait_gather(tag0[mu]);
      wait_gather(tag1[mu]);
      
      FORALLSITES(i, s) {
      
       mult_nn(&(sitez_src[i]), (matrix *)(local_pt[0][mu][i]), &tmat ); // tmat = z(n)*chi_mu(n+mu)
       scalar_mult_matrix(&tmat, s->bc[mu], &tmat);                      // tmat = bc[mu]z(n)*chi_mu(n+mu)
  
       scalar_mult_dif_matrix(&tmat,1.0,&(UpsiU[mu][i]));                // UpsiU[mu][i]=upsi -1.0*bc[mu]z(n)*chi_mu(n+mu)
        
       mult_nn_sum((matrix *)(local_pt[0][mu][i]),(matrix *)(local_pt[1][mu][i]),&(UpsiU[mu][i]) );// (UpsiU[mu][i]) = upsi -bc[mu]z(n)*chi_mu(n+mu) + chi_mu(n+mu)*Z(n+mu)
       
    
       
       scalar_mult_matrix(&(UpsiU[mu][i]) , SIGN, &(UpsiU[mu][i]) );     // (UpsiU[mu][i])  = SIGN*[upsi -bc[mu]z(n)*chi_mu(n+mu) + chi_mu(n+mu)*Z(n+mu)]
       
       scalar_mult_sum_adj_matrix(&(UpsiU[mu][i]), -1.0, &(s->f_U[mu]));
        
      }
      
      
      cleanup_gather(tag0[mu]);
      cleanup_gather(tag1[mu]);
     
    
  }
 
 
    
 #else

#endif


}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Update the momenta with the fermion force
// Assume that the multiCG has been run, with the solution in sol[j]
// Accumulate f_U for each pole into fullforce, add to momenta
// Use fullforce-->Fmunu and tempTF for temporary storage
// (Calls assemble_fermion_force, which uses many more temporaries)
double fermion_force(Real eps, Twist_Fermion *src, Twist_Fermion **sol) {
  register int i;
  register site *s;
  int mu, n;
  double returnit = 0.0;
  matrix **fullforce = malloc(sizeof(matrix*) * NUMLINK);
  matrix **varphifullforce = malloc(sizeof(matrix*) * 1);
  matrix **phifullforce = malloc(sizeof(matrix*) * 1);
  FORALLDIR(mu)
    fullforce[mu] = Fmunu[mu];    // Use Fmunu for temporary storage
    
  varphifullforce[0] = varphi_commut;
  phifullforce[0] = phi_commut ;

  // Initialize fullforce[mu]
  fermion_op(sol[0], tempTF, PLUS);
  FORALLSITES(i, s)
    scalar_mult_TF(&(tempTF[i]), amp4[0], &(tempTF[i]));
  assemble_fermion_force(sol[0], tempTF);
  FORALLDIR(mu) {
    FORALLSITES(i, s)
      adjoint(&(s->f_U[mu]), &(fullforce[mu][i]));
  }
    
    FORALLSITES(i, s){
      adjoint(&(s->f_varphi), &(varphifullforce[0][i])); //edit 
      adjoint(&(s->f_phi), &(phifullforce[0][i]));
      }
  
  for (n = 1; n < Norder; n++) {
    fermion_op(sol[n], tempTF, PLUS);
    // Makes sense to multiply here by amp4[n]...
    FORALLSITES(i, s)
      scalar_mult_TF(&(tempTF[i]), amp4[n], &(tempTF[i]));

    assemble_fermion_force(sol[n], tempTF);
    // Take adjoint but don't negate yet...
    FORALLDIR(mu) {
      FORALLSITES(i, s)
        sum_adj_matrix(&(s->f_U[mu]), &(fullforce[mu][i]));
    }
    
    
    FORALLSITES(i, s){
       sum_adj_matrix(&(s->f_varphi), &(varphifullforce[0][i])); //edit
       sum_adj_matrix(&(s->f_phi), &(phifullforce[0][i]));

        }
        
  }

  // Update the momentum from the fermion force -- sum or eps
  // Opposite sign as to gauge force,
  // because dS_G / dU = 2F_g while ds_F / dU = -2F_f
  // Move negation here as well, though adjoint remains above
  
  FORALLSITES(i, s) {
    FORALLDIR(mu) {
      scalar_mult_dif_matrix(&(fullforce[mu][i]), eps, &(s->mom[mu]));
      returnit += realtrace(&(fullforce[mu][i]), &(fullforce[mu][i]));
    }
  }
 
  
  FORALLSITES(i, s) {
    scalar_mult_dif_matrix(&(varphifullforce[0][i]), eps, &(s->mom_varphi));
    returnit += realtrace(&(varphifullforce[0][i]), &(varphifullforce[0][i]));
  }

  FORALLSITES(i, s) {
    scalar_mult_dif_matrix(&(phifullforce[0][i]), eps, &(s->mom_phi));
    returnit += realtrace(&(phifullforce[0][i]), &(phifullforce[0][i]));
  }
  

  
  g_doublesum(&returnit);

  free(fullforce);
  free(varphifullforce);
  free(phifullforce);

  // Reset Fmunu, varphi_commutator
  compute_Fmunu();
  compute_phi_commutator();
  compute_varphi_commutator();
  return (eps * sqrt(returnit) / volume);
}
#endif
// -----------------------------------------------------------------
