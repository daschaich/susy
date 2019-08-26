// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Dot product of non-diagonal gamma matrices with scalars
// Used by Yukawa term in fermion operator
// (The last two gamma matrices are diagonal in our basis)
void build_Gamma_X() {
  register int i, j, k, l;
  register site *s;
  Real tr;

  for (j = 0; j < NCHIRAL_FERMION; j++) {
    for (k = 0; k < NCHIRAL_FERMION; k++) {
      // Overwrite Gamma_X[j][k][i], potentially with zero
      tr = Gamma[0][j][k];
      FORALLSITES(i, s)
        scalar_mult_matrix(&(s->X[0]), tr, &(Gamma_X[j][k][i]));

      // Add remaining non-zero contributions
      for (l = 1; l < NSCALAR - 2; l++) {
        tr = Gamma[l][j][k];
        if (tr * tr > SQ_TOL) {
          FORALLSITES(i, s)
            scalar_mult_sum_matrix(&(s->X[l]), tr, &(Gamma_X[j][k][i]));
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void Dplus(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  register int i, j, k;
  register site *s;
  matrix tmat;
  msg_tag *tag[NCHIRAL_FERMION];

  for (j = 0; j < NCHIRAL_FERMION; j++) {
    k = j + NCHIRAL_FERMION;
    tag[j] = start_gather_field(src[k], sizeof(matrix),
                                TUP, EVENANDODD, gen_pt[j]);
  }

  for (j = 0; j < NCHIRAL_FERMION; j++) {     // Overwrite dest
    k = j + NCHIRAL_FERMION;
    wait_gather(tag[j]);
    FORALLSITES(i, s) {
      mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
      scalar_mult_na(&tmat, &(s->link), s->bc[0], &(dest[j][i]));
      dif_matrix(&(src[k][i]), &(dest[j][i]));
    }
    cleanup_gather(tag[j]);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use tempmat for temporary storage
void Dminus(matrix *src[NFERMION], matrix *dest[NFERMION]) {
  register int i, j, k;
  register site *s;
  matrix tmat;
  msg_tag *tag;

  // Include negative sign in product that is then gathered
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    k = j + NCHIRAL_FERMION;
    FORALLSITES(i, s) {
      mult_an(&(s->link), &(src[j][i]), &tmat);
      scalar_mult_nn(&tmat, &(s->link), -1.0, &(tempmat[i]));
    }
    tag = start_gather_field(tempmat, sizeof(matrix),
                             TDOWN, EVENANDODD, gen_pt[0]);

    wait_gather(tag);
    FORALLSITES(i, s) {           // Overwrite dest
      scalar_mult_add_matrix(&(src[j][i]), (matrix *)(gen_pt[0][i]),
                             s->bc[1], &(dest[k][i]));
    }
    cleanup_gather(tag);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
// Uses tempmat for temporary storage through Dminus
#ifndef PUREGAUGE
void fermion_op(matrix *src[NFERMION], matrix *dest[NFERMION], int sign) {
  register int i, j, k, m, n;
  register site *s;
  Real tr;
  // Sqrt factor from Eq. 14 in 2 Jun 2019 notes
  complex tc = cmplx(0.0, 1.0 / sqrt(2.0));
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
    m = j + NCHIRAL_FERMION;
    for (k = 0; k < NCHIRAL_FERMION; k++) {
      n = k + NCHIRAL_FERMION;
      tr = sign * mass_fermion * Gamma123[j][k];      // sign = +/- 1
      if (tr * tr > SQ_TOL) {
        FORALLSITES(i, s) {
          scalar_mult_sum_matrix(&(src[n][i]), tr, &(dest[j][i]));
          scalar_mult_dif_matrix(&(src[k][i]), tr, &(dest[m][i]));
        }
      }
    }
  }
#endif

  // Yukawa terms -- assume build_Gamma_X has already been run
  // Sqrt factor from Eq. 14 in 2 Jun 2019 notes
  // Gamma_X[j][k] embedded into -i*sigma_2 = ( 0 -1 \\ 1 0 ) in 8x8 blocks
  tr = (Real)sign / sqrt(2.0);
  for (j = 0; j < NCHIRAL_FERMION; j++) {
    m = j + NCHIRAL_FERMION;
    for (k = 0; k < NCHIRAL_FERMION; k++) {
      n = k + NCHIRAL_FERMION;

      // Can skip j==k since Gamma_X[j][k] anti-symmetric
      if (j == k)
        continue;

      FORALLSITES(i, s) {
        // Upper-right block takes n to j with negative sign
        mult_nn(&(Gamma_X[j][k][i]), &(src[n][i]), &tmat);
        mult_nn_dif(&(src[n][i]), &(Gamma_X[j][k][i]), &tmat);
        scalar_mult_dif_matrix(&tmat, tr, &(dest[j][i]));

        // Lower-left block takes k to m
        mult_nn(&(Gamma_X[j][k][i]), &(src[k][i]), &tmat);
        mult_nn_dif(&(src[k][i]), &(Gamma_X[j][k][i]), &tmat);
        scalar_mult_sum_matrix(&tmat, tr, &(dest[m][i]));
      }
    }

    // Last 2 gammas are diagonal
    k = NSCALAR - 2;
    n = NSCALAR - 1;
    FORALLSITES(i, s) {
      // Gamma[7] = ( -1 0 \\ 0 1 ) in 8x8 blocks
      mult_nn(&(s->X[k]), &(src[j][i]), &tmat);
      mult_nn_dif(&(src[j][i]), &(s->X[k]), &tmat);
      scalar_mult_dif_matrix(&tmat, tr, &(dest[j][i]));

      mult_nn(&(s->X[k]), &(src[m][i]), &tmat);
      mult_nn_dif(&(src[m][i]), &(s->X[k]), &tmat);
      scalar_mult_sum_matrix(&tmat, tr, &(dest[m][i]));

      // Gamma[8] = ( -i 0 \\ 0 -i ) in 8x8 blocks
      // tc purely imaginary --> same for both sign = +/-1
      mult_nn(&(s->X[n]), &(src[j][i]), &tmat);
      mult_nn_dif(&(src[j][i]), &(s->X[n]), &tmat);
      c_scalar_mult_dif_mat(&tmat, &tc, &(dest[j][i]));

      mult_nn(&(s->X[n]), &(src[m][i]), &tmat);
      mult_nn_dif(&(src[m][i]), &(s->X[n]), &tmat);
      c_scalar_mult_dif_mat(&tmat, &tc, &(dest[m][i]));
    }
  }

  // Overall factor of 1/2
  tr = 0.5 * kappa;
  for (j = 0; j < NFERMION; j++) {
    FORALLSITES(i, s)
      scalar_mult_matrix(&(dest[j][i]), tr, &(dest[j][i]));
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
