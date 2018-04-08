// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "susy_includes.h"
// -----------------------------------------------------------------

void build_Gamma_X() {
  register int i,j,k,l;
  register site *s;
  
  for(j=0; j<NCHIRAL_FERMION; j++) {
    for(k=0; k<NCHIRAL_FERMION; k++) {
      FORALLSITES(i, s) {
        clear_mat(&(Gamma_X[j][k][i]));
        for(l=0; l<NSCALAR-2; l++) {
          scalar_mult_sum_matrix(&(s->X[l]), Gamma[l].e[j][k],
                                 &(Gamma_X[j][k][i]));
        }
      }
    }
  }
}

void Dplus(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  register int i,j;
  register site *s;
  matrix tmat;
  msg_tag *tag[NCHIRAL_FERMION];
  
  for(j=0; j<NCHIRAL_FERMION; j++) {
    tag[j] = start_gather_field(src[j+NCHIRAL_FERMION], sizeof(matrix),
                                TUP, EVENANDODD, gen_pt[j]);
  }
  
  for(j=0; j<NCHIRAL_FERMION; j++) {
    wait_gather(tag[j]);
    // Initialize dest[i]
    FORALLSITES(i, s) {
      mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
      scalar_mult_na(&tmat, &(s->link), s->bc[0], &(dest[j][i]));
      dif_matrix(&(src[j+NCHIRAL_FERMION][i]), &(dest[j][i]));
    }
    cleanup_gather(tag[j]);
  }
}

// This uses tempmat for temporary storage
void Dminus(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  register int i,j;
  register site *s;
  matrix tmat;
  msg_tag *tag;
  
  for(j=0; j<NCHIRAL_FERMION; j++) {
    // Initialize dest[i]
    FORALLSITES(i, s) {
      mult_an(&(s->link), &(src[j][i]), &tmat);
      mult_nn(&tmat, &(s->link), &(tempmat[i]));
    }

    tag = start_gather_field(tempmat, sizeof(matrix),
                                TDOWN, EVENANDODD, gen_pt[0]);
    wait_gather(tag);
    // Initialize dest[i]
    FORALLSITES(i, s) {
      scalar_mult_add_matrix(&(src[j][i]), (matrix *)(gen_pt[0][i]),
                             -1.0*s->bc[1], &(dest[j+NCHIRAL_FERMION][i]));
    }
    cleanup_gather(tag);
  }
}


// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
#ifndef PUREGAUGE
void fermion_op(matrix *src[NFERMION], matrix *dest[NFERMION], int sign) {
  register int i,j,k;
  register site *s;
  Real tr;
  complex tc = cmplx(0.0, 1.0);
  matrix tmat;

  // Fermion kinetic term
  Dplus(src, dest);
  Dminus(src, dest);
  
  if (sign == -1) {
    for (j = 0; j < NFERMION; j++) {
      FORALLSITES(i, s)
        scalar_mult_matrix(&(dest[j][i]), -1.0, &(dest[j][i]));
    }
  }
  
#ifdef BMN
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    for (k = 0; k < NCHIRAL_FERMION; k++) {
      tr = 0.75 * sign * mu * Gamma123.e[j][k]; // sign = +/- 1
      FORALLSITES(i, s) {
        scalar_mult_sum_matrix(&(src[k+NCHIRAL_FERMION][i]), tr, &(dest[j][i]));
        scalar_mult_dif_matrix(&(src[k][i]), tr, &(dest[j+NCHIRAL_FERMION][i]));
      }
    }
  }
#endif
  
  // Assume build_Gamma_X has already been run
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    for (k = 0; k < NCHIRAL_FERMION; k++) {
      FORALLSITES(i, s) {
        mult_nn(&(Gamma_X[j][k][i]), &(src[k+NCHIRAL_FERMION][i]), &tmat);
        mult_nn_dif(&(src[k+NCHIRAL_FERMION][i]), &(Gamma_X[j][k][i]), &tmat);
        scalar_mult_sum_matrix(&tmat, (Real)sign, &(dest[j][i]));
        
        mult_nn(&(Gamma_X[k][j][i]), &(src[k][i]), &tmat);
        mult_nn_dif(&(src[k][i]), &(Gamma_X[k][j][i]), &tmat);
        scalar_mult_sum_matrix(&tmat, (Real)sign, &(dest[j+NCHIRAL_FERMION][i]));
      }
    }
    // Last 2 gammas are diagonal
    FORALLSITES(i, s) {
      mult_nn(&(s->X[NSCALAR-2]), &(src[j][i]), &tmat);
      mult_nn_dif(&(src[j][i]), &(s->X[NSCALAR-2]), &tmat);
      scalar_mult_dif_matrix(&tmat, (Real)sign, &(dest[j][i]));
      
      mult_nn(&(s->X[NSCALAR-2]), &(src[j+NCHIRAL_FERMION][i]), &tmat);
      mult_nn_dif(&(src[j+NCHIRAL_FERMION][i]), &(s->X[NSCALAR-2]), &tmat);
      scalar_mult_sum_matrix(&tmat, (Real)sign, &(dest[j+NCHIRAL_FERMION][i]));
      
      mult_nn(&(s->X[NSCALAR-1]), &(src[j][i]), &tmat);
      mult_nn_dif(&(src[j][i]), &(s->X[NSCALAR-1]), &tmat);
      c_scalar_mult_sum_mat(&tmat, &tc, &(dest[j][i]));
      
      mult_nn(&(s->X[NSCALAR-1]), &(src[j+NCHIRAL_FERMION][i]), &tmat);
      mult_nn_dif(&(src[j+NCHIRAL_FERMION][i]), &(s->X[NSCALAR-1]), &tmat);
      c_scalar_mult_sum_mat(&tmat, &tc, &(dest[j+NCHIRAL_FERMION][i]));
      
      scalar_mult_matrix(&(dest[j][i]), 0.5, &(dest[j][i]));
      scalar_mult_matrix(&(dest[j+NCHIRAL_FERMION][i]), 0.5,
                         &(dest[j+NCHIRAL_FERMION][i]));
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Squared matrix--vector operation
//   dest = (D^2).src
// Use temp_ferm for temporary storage
#ifndef PUREGAUGE
void DSq(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  fermion_op(src, temp_ferm, PLUS);
  fermion_op(temp_ferm, dest, MINUS);
}
#endif
// -----------------------------------------------------------------
