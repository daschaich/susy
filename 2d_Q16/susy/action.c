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
  complex tc, tc1;
  msg_tag *mtag0 = NULL;

  FORALLDIR(mu) {
    FORALLSITES(i, s) {
      mult_na(&(s->link[mu]), &(s->link[mu]), &(tempmat[i])); // tempmat[i]= U_mu(i)U_mu_bar(i)
      mult_an(&(s->link[mu]), &(s->link[mu]), &(tempmat2[i])); // tempmat2[i]= U_mu_bar(i)U_mu(i)
    }

    // Gather tempmat2 from below
    mtag0 = start_gather_field(tempmat2, sizeof(matrix),
                               goffset[mu] + 1, EVENANDODD, gen_pt[0]); //gen_pt[0]= U_mu_bar(i-mu)U_mu(i-mu)
    wait_gather(mtag0);
    if (mu == 0) {
      FORALLSITES(i, s)        // Initialize
        sub_matrix(&(tempmat[i]), (matrix *)(gen_pt[0][i]), &(DmuUmu[i])); //DmuUmu(i)= U_0(i)U_mu_bar(i) - U_0_bar(i-0)U_0(i-0)
    }
    else {
      FORALLSITES(i, s) {
        sum_matrix(&(tempmat[i]), &(DmuUmu[i]));   // DmuUmu[i] = DmuUmu[i] + U_mu(i)U_mu_bar(i) = U_0(i)U_mu_bar(i) - U_0_bar(i-0)U_0(i-0) + U_mu(i)U_mu_bar(i)
        dif_matrix((matrix *)(gen_pt[0][i]), &(DmuUmu[i])); // DmuUmu[i] = DmuUmu[i] - gen_pt[0][i] = U_0(i)U_mu_bar(i) - U_0_bar(i-0)U_0(i-0) + U_mu(i)U_mu_bar(i) - U_mu_bar(i-mu)U_mu(i-mu)
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
  int mu, nu, index;
  msg_tag *mtag0 = NULL, *mtag1 = NULL;
  complex tc;

  FORALLDIR(mu) {
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      index = plaq_index[mu][nu];

      mtag0 = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                goffset[mu], EVENANDODD, gen_pt[0]); //gen_pt[0]=U_nu(n+mu)
      mtag1 = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                                goffset[nu], EVENANDODD, gen_pt[1]); //gen_pt[1]=U_mu(n+nu)
      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES(i, s) {
        mult_nn(&(s->link[mu]), (matrix *)(gen_pt[0][i]), &(tempmat[i])); // U_mu(n)U_nu(n+mu)
        mult_nn(&(s->link[nu]), (matrix *)(gen_pt[1][i]), &(tempmat2[i])); //U_nu(n) U_mu(n+nu)
        sub_matrix(&(tempmat[i]), &(tempmat2[i]), &(Fmunu[index][i]));  // U_mu(n)U_nu(n+mu) - U_nu(n) U_mu(n+nu)
        
      
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
}
// -----------------------------------------------------------------


// D+ original
void compute_Dmuphi() {
  register int i, mu;
  register site *s;
  msg_tag *tag[NUMLINK],*tag0[0]; // , *mtag1 = NULL  ;

  FORALLDIR(mu) {
  
   tag0[0] = start_gather_site(F_OFFSET(phi), sizeof(matrix),
                              goffset[mu], EVENANDODD, gen_pt[0]);// gen_pt[0] = phi(n+mu)


    wait_gather(tag0[0]);
      FORALLSITES(i, s) {
      
       mult_nn(&(s->link[mu]), (matrix *)(gen_pt[0][i]), &(tempmat[i]));  // c <-- a * b => tempmat[i] = U_mu(i)phi(i+mu)
       mult_nn(&(s->phi), &(s->link[mu]),&(tempmat2[i]));                  // c <-- a * b => tempmat2[i] = phi(i)U_mu(i)
       sub_matrix(&(tempmat[i]), &(tempmat2[i]), &(Dmuphi[mu][i]));             // c <-- a - b => Dmuphi[mu][i] = U_mu(n)phi(n+mu) - phi(n)U_mu(n)
                                                                            // not c <-- c - a - b  
       }
      cleanup_gather(tag0[0]);
     
    
  }
}




//D+ original 

void compute_Dmuvarphi() {
  register int i, mu;
  register site *s;
  msg_tag *tag[NUMLINK],*tag0[0]; // , *mtag1 = NULL  ;

  FORALLDIR(mu) {
  
   tag0[0] = start_gather_site(F_OFFSET(varphi), sizeof(matrix),
                              goffset[mu], EVENANDODD, gen_pt[0]);// gen_pt[0] = varphi(n+mu)


    wait_gather(tag0[0]);
      FORALLSITES(i, s) {
      
       mult_nn(&(s->link[mu]), (matrix *)(gen_pt[0][i]), &(tempmat[i]));  // c <-- a * b => tempmat[i] = U_mu(i)phi(i+mu)
       mult_nn(&(s->varphi), &(s->link[mu]),&(tempmat2[i]));                  // c <-- a * b => tempmat2[i] = phi(i)U_mu(i)
       sub_matrix(&(tempmat[i]), &(tempmat2[i]), &(Dmuvarphi[mu][i]));             // c <-- a - b => Dmuphi[mu][i] = U_mu(n)phi(n+mu) - phi(n)U_mu(n)
                                                                            // not c <-- c - a - b  
       }
      cleanup_gather(tag0[0]);
     
    
  }
}



//----[varphi,phi]

void compute_varphi_phi_commutator() 
{
  register int i;
  register site *s;

 // msg_tag *mtag0 = NULL, *mtag1 = NULL;


      FORALLSITES(i, s) 
      {
      
        
       
        mult_nn(&(s->varphi),&(s->phi), &(tempmat[i])) ;    // c <-- a * b = varphi(n)phi(n)
        
        mult_nn(&(s->phi), &(s->varphi), &(tempmat2[i]));  // c <-- a * b = phi(n)varphi(n)
        
        sub_matrix(&(tempmat[i]), &(tempmat2[i]),&(vp_commut[i]));  // c <-- a - b = phi(n)varphi(n) - varphi(n)phi(n)
       
      }

    
  
}


//----[phi^bar,phi]

void compute_phi_commutator() 
{
  register int i;
  register site *s;
  complex tc,tc1,tc2;
  matrix tmat , tmat2;

 // msg_tag *mtag0 = NULL, *mtag1 = NULL;


      FORALLSITES(i, s) 
      {
      
        //clear_mat(&(phi_commut[i]));
      
        mult_an(&(s->phi), &(s->phi), &(tempmat[i]));  // c <-- adag * b  = phi^bar(n)phi(n)
  
        mult_na(&(s->phi), &(s->phi), &(tempmat2[i]));  // c <-- a * bdag  = phi(n)phi^bar(n)
        
      
       sub_matrix(&(tempmat[i]), &(tempmat2[i]),&(phi_commut[i]));  // c <-- a - b = phi^bar(n)phi(n) - phi(n)phi^bar(n)
        
      
      }

      
}



//----[varphi^bar,varphi]

void compute_varphi_commutator() 
{
  register int i;
  register site *s;

  //msg_tag *mtag0 = NULL, *mtag1 = NULL;


      FORALLSITES(i, s) 
      {
      

      
        mult_an(&(s->varphi), &(s->varphi), &(tempmat[i]));  // c <-- adag * b  = varphi^bar(n)varphi(n)
        
        mult_na(&(s->varphi), &(s->varphi), &(tempmat2[i]));  // c <-- a * bdag  = varphi(n)varphi^bar(n)
       
        sub_matrix(&(tempmat[i]), &(tempmat2[i]),&(varphi_commut[i]));  // c <-- a - b = varphi^bar(n)varphi(n) - varphi(n)varphi^bar(n)
       
      }

    
  
}




//-----------------------------end 7/11/23--------------------------





// -----------------------------------------------------------------
// Standard gauge contribution to the action
// Include tunable coefficient C2 in the d^2 term
double gauge_action(int do_det) {
  register int i;
  register site *s;
  int index;
  double g_action = 0.0, norm = 0.5 * C2;
  complex tc;
  matrix tmat, tmat2;

  FORALLSITES(i, s) {
    // d^2 term normalized by C2 / 2
    // DmuUmu includes the plaquette determinant contribution if G is non-zero
    mult_nn(&(DmuUmu[i]), &(DmuUmu[i]), &tmat);
    scalar_mult_matrix(&tmat, norm, &tmat);
    

    // F^2 term
    
    
    
    for (index = 0; index < NPLAQ; index++) {
      mult_an(&(Fmunu[index][i]), &(Fmunu[index][i]), &tmat2);
      scalar_mult_sum_matrix(&tmat2, 2.0, &tmat); //2
    }
        
    
    //---------edited 7/9/23 end----------------

    if (do_det == 1) {
      det_project(&tmat, &tmat2);
      tc = trace(&tmat2);
    }
    else
      tc = trace(&tmat);

    g_action += tc.real;
    
   // printf("node = %d site = %d trace gauge = %f\n",this_node,i,tc.real); //FORALLSITE(i,s)=for(i=0,s=lattice;i<sites_on_node;i++,s++)
    
   // printf("trace gauge = %d\n",tc.real);
#ifdef DEBUG_CHECK
    if (fabs(tc.imag) > IMAG_TOL)
      printf("node%d WARNING: Im[s_B[%d]] = %.4g\n", this_node, i, tc.imag);
#endif
  }
  g_action *= kappa;      // kappa = N/4*lambda
  g_doublesum(&g_action); // Sum a double over all nodes of mpi
  return g_action;
}
// -----------------------------------------------------------------

//---edited-----------------
// scalar action 



double sa(int do_det) {
  register int i;
  register site *s;
  double s_a = 0.0;
  int index;
  complex tc;
  matrix tmat, tmat2;
  

  FORALLSITES(i, s) {
    
    clear_mat(&tmat);
   
   
   
    
  // 1/2[phi^bar,phi]^2 term 
    
  //  mult_nn(&(phi_commut[i]), &(phi_commut[i]), &tmat2); // c <-- a * b
  //  scalar_mult_sum_matrix(&tmat2,0.5, &tmat); // c <-- c + s * b
  
  
   
   mult_nn(&(phi_commut[i]), &(phi_commut[i]), &tmat); // c <-- a * b
   scalar_mult_matrix(&tmat, 0.5, &tmat);
   
   
    // 1/2[varphi^bar,varphi]^2 term 
    
    
    // mult_nn(&(varphi_commut[i]), &(phi_commut[i]), &tmat); // c <-- a * b
    // scalar_mult_matrix(&tmat, 0.5, &tmat);
    
    
   
    mult_nn(&(varphi_commut[i]), &(varphi_commut[i]), &tmat2); // c <-- a * b
    scalar_mult_sum_matrix(&tmat2, 0.5, &tmat);
    
    
    
    // DmuUmu*[phi^bar,phi] term --->FD + Incorrect force
    
    mult_nn(&(DmuUmu[i]),&(phi_commut[i]), &tmat2); // c <-- a * b
    scalar_mult_sum_matrix(&tmat2, 1.0, &tmat); // c <-- c + s * b
 
       
    
   // DmuUmu*[varphi^bar,varphi] term --->FD + incorrect force
    
    mult_nn(&(DmuUmu[i]),&(varphi_commut[i]), &tmat2); // c <-- a * b
    scalar_mult_sum_matrix(&tmat2, 1.0, &tmat); // c <-- c + s * b
    
    
    
    // -2|Dmuphi|^2 term  --->FD + Incorrect force
    
    for (index = 0; index < NUMLINK; index++) {
    
    mult_an(&(Dmuphi[index][i]), &(Dmuphi[index][i]), &tmat2); // c <-- adag * b => tmat2 = Dmuphi^bar[index][i]*Dmuphi[index][i]
    scalar_mult_sum_matrix(&tmat2, 2.0, &tmat); // c <-- c + s * b  => tmat = tmat + Dmuphi^bar[index][i]*Dmuphi[index][i]
  
    }
    
    // -2|Dmuvarphi|^2 term  --->FD+Incorrect force
    
    for (index = 0; index < NUMLINK; index++) {
    
    mult_an(&(Dmuvarphi[index][i]), &(Dmuvarphi[index][i]), &tmat2); // c <-- adag * b => tmat2 = Dmuvarphi^bar[index][i]*Dmuvarphi[index][i]
    scalar_mult_sum_matrix(&tmat2, 2.0, &tmat); // c <-- c + s * b  => tmat = tmat + Dmuvarphi^bar[index][i]*Dmuvarphi[index][i]
  
    }
    

    
    // [phi^bar,phi]*[varphi^bar,varphi] term FD + incorrect force 
    
      
    mult_nn(&(phi_commut[i]), &(varphi_commut[i]), &tmat2); // c <-- a * b
    scalar_mult_sum_matrix(&tmat2, 1.0, &tmat);
    
    
   
    
    
    // -2|[varphi,phi]|^2 term-----> FD + incorrect Force
    
    mult_an(&(vp_commut[i]), &(vp_commut[i]), &tmat2); // c <-- adag * b => tmat2 = Dmuvarphi^bar[index][i]*Dmuvarphi[index][i]
    scalar_mult_sum_matrix(&tmat2, 2.0, &tmat); // c <-- c + s * b  => tmat = tmat + Dmuvarphi^bar[index][i]*Dmuvarphi[index][i]
    
    
    
    //---------7/9/23 end----------------
    
    
    if (do_det == 1) {
      det_project(&tmat, &tmat2);
      tc = trace(&tmat2);
    }

    else
      tc = trace(&tmat);

    s_a += tc.real;
    
  //  printf("node = %d site = %d trace scalar = %f\n",this_node,i,tc.real);
    
#ifdef DEBUG_CHECK
    if (fabs(tc.imag) > IMAG_TOL)
      printf("node%d WARNING: Im[s_B[%d]] = %.4g\n", this_node, i, tc.imag);
#endif
  }
  s_a *= kappa;      // kappa = N/4*lambda
  g_doublesum(&s_a); // Sum a double over all nodes of mpi
  return s_a;
}



//term action 




// flat direction


double bmass_action() {
  register int i, a;
  register site *s;
  double sum = 0.0;
#ifdef EIG_POT
  matrix tmat;
#else
  Real tr;
#endif

  FORALLSITES(i, s)   {
  
      //FORALLDIR(a) {
#ifdef EIG_POT

      //Tr [U_a(x) Ubar_a(x) - I]^2
      
      
      FORALLDIR(a) {
      
      mult_na(&(s->link[a]), &(s->link[a]), &tmat);
      scalar_add_diag(&tmat, -1.0);
      sum += realtrace(&tmat, &tmat); // tr([U_a(x) Ubar_a(x) - I]^dag*[U_a(x) Ubar_a(x) - I])
      
                   }
       
      
      // Tr[phi(x)phi^bar(x)]^2    
      mult_na(&(s->phi), &(s->phi), &tmat);   
      scalar_add_diag(&tmat, -1.0);      
      sum += realtrace(&tmat, &tmat);       //**correct*** not Tr[tmat^bar*tmat] but   Tr[tmat]    
      
               
      // Tr[varphi(x)varphi^bar(x)]^2    
      mult_na(&(s->varphi), &(s->varphi), &tmat);   
      scalar_add_diag(&tmat, -1.0);      
      sum += realtrace(&tmat, &tmat);       //**correct*** not Tr[tmat^bar*tmat] but   Tr[tmat]        
      
      
                  
      
#else

      FORALLDIR(a) {
      // tr = one_ov_N * realtrace(&(s->link[a]), &(s->link[a])) + one_ov_N * realtrace(&(s->phi), &(s->phi)) + one_ov_N * realtrace(&(s->varphi), &(s->varphi)) - 1.0;
      tr = one_ov_N * realtrace(&(s->link[a]), &(s->link[a])) - 1.0;
      sum += tr * tr;
                   }
                   
     tr = one_ov_N * realtrace(&(s->phi), &(s->phi)) - 1.0 ;
     sum += tr * tr;
     tr = one_ov_N * realtrace(&(s->varphi), &(s->varphi)) - 1.0 ;
     sum += tr * tr;             


#endif              

                     }
                   
  sum *= kappa * bmass * bmass;
  g_doublesum(&sum);
  return sum;
}

 



// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette determinant contribution to the action
// Assume compute_plaqdet() has already been run


//------- to be edited under DR -------
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
  det_action *= kappa_u1; // Plaquette determinant coupling 'k'
  g_doublesum(&det_action);
  return det_action;
}

//-----end-------------


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
      sum += (double)realtrace(&(s->mom[mu]), &(s->mom[mu])); // trace (adag*b)
     
    sum += (double)realtrace(&(s->mom_phi), &(s->mom_phi));
    sum += (double)realtrace(&(s->mom_varphi), &(s->mom_varphi));
  }
  //edited---------
  /*
    FORALLSITES(i, s) {
      sum += (double)realtrace(&(s->mom_phi), &(s->mom_phi));
        }
  */
  
  
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Print out zeros for pieces of the action that aren't included
double action(Twist_Fermion **src, Twist_Fermion ***sol) {
  double g_act, sca, BA, bmass_act = 0.0,sp_act, p_act, f_act, det_act = 0.0;
  double total;
  register int m;
  register site *s;
  g_act = gauge_action(NODET);
  sca   = sa(NODET);
 // BA    = g_act+sca;

  
  if (bmass > IMAG_TOL)
    bmass_act = bmass_action();
  if (kappa_u1 > IMAG_TOL)
    det_act = det_action();

  node0_printf("action: ");
  node0_printf("gauge %.8g bmass %.8g scalar %.8g ", g_act, bmass_act, sca);
 
  //node0_printf("gauge %.8g bmass %.8g ", g_act, bmass_act); 
  node0_printf("det %.8g ", det_act);

  total = g_act + bmass_act + det_act + sca;
 
 //total = g_act + bmass_act + det_act;
  //total =  sca;
  
 //node0_printf("scalar_action = %.8g ", sca);
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
