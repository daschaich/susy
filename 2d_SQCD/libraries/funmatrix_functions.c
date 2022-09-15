// A file that containts the functions needed for funmatrix operations.
// TODO split in different files at the end.
// TODO add all functions in include/susy.h

//Currently funmatrices are NCOL x NCOLF matrices
//and funamatrices are NCOLF x NCOL matrices

//For all
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

//For fun_dumpmat and funa_dumpmat
#include <stdio.h>

#define DEBUGFUNMATRIX

// -----------------------------------------------------------------
// Add funmatrices
// c <-- c + b
// c <-- a + b

void fun_sum_matrix(funmatrix *b, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real += b->e[i][j].real;
      c->e[i][j].imag += b->e[i][j].imag;
    }
  }
}

void fun_add_matrix(funmatrix *a, funmatrix *b, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + b->e[i][j].imag;
    }
  }
}

// Add funamatrices
// c <-- c + b
// c <-- a + b

void funa_sum_matrix(funamatrix *b, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[i][j].real;
      c->e[i][j].imag += b->e[i][j].imag;
    }
  }
}

void funa_add_matrix(funamatrix *a, funamatrix *b, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Adjoint of a funmatrix
// b <-- (+/-)adag

void fun_adjoint(funmatrix *a, funamatrix *b) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++)
      CONJG(a->e[j][i], b->e[i][j]);
  }
}

void fun_neg_adjoint(funmatrix *a, funamatrix *b) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = -a->e[j][i].real;
      b->e[i][j].imag = a->e[j][i].imag;
    }
  }
}

// Adjoint of a funamatrix

void funa_adjoint(funamatrix *a, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++)
      CONJG(a->e[j][i], b->e[i][j]);
  }
}

void funa_neg_adjoint(funamatrix *a, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      b->e[i][j].real = -a->e[j][i].real;
      b->e[i][j].imag = a->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Clear the given funmatrix

void fun_clear_mat(funmatrix *m) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      m->e[i][j].real = 0.0;
      m->e[i][j].imag = 0.0;
    }
  }
}

// Clear the given funamatrix

void funa_clear_mat(funamatrix *m) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      m->e[i][j].real = 0.0;
      m->e[i][j].imag = 0.0;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Add or subtract result of adjoint of complex scalar multiplication on funmatrix
// c <-- c + (s * b)dag
// c <-- c - (s * b)dag
// (s * b).real = s.real * b.real - s.imag * b.imag
// (s * b).imag = s.imag * b.real + s.real * b.imag
// Then flip i<-->j and negate (s * b).imag

void fun_c_scalar_mult_sum_adj_mat(funmatrix *b, complex *s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[j][i].real * s->real - b->e[j][i].imag * s->imag;
      c->e[i][j].imag -= b->e[j][i].imag * s->real + b->e[j][i].real * s->imag;
    }
  }
}

void fun_c_scalar_mult_dif_adj_mat(funmatrix *b, complex *s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[j][i].real * s->real - b->e[j][i].imag * s->imag;
      c->e[i][j].imag += b->e[j][i].imag * s->real + b->e[j][i].real * s->imag;
    }
  }
}

// Add or subtract result of adjoint of complex scalar multiplication on funamatrix

void funa_c_scalar_mult_sum_adj_mat(funamatrix *b, complex *s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real += b->e[j][i].real * s->real - b->e[j][i].imag * s->imag;
      c->e[i][j].imag -= b->e[j][i].imag * s->real + b->e[j][i].real * s->imag;
    }
  }
}

void funa_c_scalar_mult_dif_adj_mat(funamatrix *b, complex *s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real -= b->e[j][i].real * s->real - b->e[j][i].imag * s->imag;
      c->e[i][j].imag += b->e[j][i].imag * s->real + b->e[j][i].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Add result of complex scalar multiplication on funmatrix
// c <-- c + s * b

void fun_c_scalar_mult_sum_mat(funmatrix *b, complex *s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real += b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag += b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}

// Add result of complex scalar multiplication on funamatrix

void funa_c_scalar_mult_sum_mat(funamatrix *b, complex *s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag += b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Add result of complex scalar multiplication on adjoint funmatrix
// c <-- c + s * bdag
// Just flip i<-->j and negate b.imag

void fun_c_scalar_mult_sum_mat_adj(funmatrix *b, complex *s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[j][i].real * s->real + b->e[j][i].imag * s->imag;
      c->e[i][j].imag -= b->e[j][i].imag * s->real - b->e[j][i].real * s->imag;
    }
  }
}

// Add result of complex scalar multiplication on adjoint funamatrix

void funa_c_scalar_mult_sum_mat_adj(funamatrix *b, complex *s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real += b->e[j][i].real * s->real + b->e[j][i].imag * s->imag;
      c->e[i][j].imag -= b->e[j][i].imag * s->real - b->e[j][i].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Complex scalar multiplication on funmatrix
// c <-- s * b

void fun_c_scalar_mult_mat(funmatrix *b, complex *s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real = b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag = b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}

// Complex scalar multiplication on funamatrix

void funa_c_scalar_mult_mat(funamatrix *b, complex *s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag = b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Subtract result of complex scalar multiplication on funmatrix
// c <-- c - s * b

void fun_c_scalar_mult_dif_mat(funmatrix *b, complex *s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real -= b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag -= b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}

// Subtract result of complex scalar multiplication on funmatrix

void funa_c_scalar_mult_dif_mat(funamatrix *b, complex *s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag -= b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Print the given funmatrix

void fun_dumpmat(funmatrix *m) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++)
      printf("  (%.4g, %.4g)", m->e[i][j].real, m->e[i][j].imag);
    printf("\n");
  }
}

// Print the given funamatrix

void funa_dumpmat(funamatrix *m) {
  int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++)
      printf("  (%.4g, %.4g)", m->e[i][j].real, m->e[i][j].imag);
    printf("\n");
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Matrix multiplication with adjoint of first matrix
// c <-- c + adag * b
// c <-- c - adag * b
// c <-- adag * b
// a,b - funamatrices

void funa_mult_an_sum(funamatrix *a, funamatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real += a->e[k][i].real * b->e[k][j].real
                         + a->e[k][i].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[k][i].real * b->e[k][j].imag
                         - a->e[k][i].imag * b->e[k][j].real;
      }
    }
  }
}

void funa_mult_an_dif(funamatrix *a, funamatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real -= a->e[k][i].real * b->e[k][j].real
                         + a->e[k][i].imag * b->e[k][j].imag;
        c->e[i][j].imag -= a->e[k][i].real * b->e[k][j].imag
                         - a->e[k][i].imag * b->e[k][j].real;
      }
    }
  }
}

void funa_mult_an(funamatrix *a, funamatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = a->e[0][i].real * b->e[0][j].real
                      + a->e[0][i].imag * b->e[0][j].imag;
      c->e[i][j].imag = a->e[0][i].real * b->e[0][j].imag
                      - a->e[0][i].imag * b->e[0][j].real;
      for (k = 1; k < NCOLF; k++) {
        c->e[i][j].real += a->e[k][i].real * b->e[k][j].real
                         + a->e[k][i].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[k][i].real * b->e[k][j].imag
                         - a->e[k][i].imag * b->e[k][j].real;
      }
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Matrix multiplication with adjoint of second matrix
// c <-- c + a * bdag
// c <-- c - a * bdag
// c <-- a * bdag
// a,b - funmatrices

void fun_mult_na_sum(funmatrix *a, funmatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[j][k].real
                         + a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[j][k].real
                         - a->e[i][k].real * b->e[j][k].imag;
      }
    }
  }
}

void fun_mult_na_dif(funmatrix *a, funmatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real -= a->e[i][k].real * b->e[j][k].real
                         + a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag -= a->e[i][k].imag * b->e[j][k].real
                         - a->e[i][k].real * b->e[j][k].imag;
      }
    }
  }
}

void fun_mult_na(funmatrix *a, funmatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = a->e[i][0].real * b->e[j][0].real
                      + a->e[i][0].imag * b->e[j][0].imag;
      c->e[i][j].imag = a->e[i][0].imag * b->e[j][0].real
                      - a->e[i][0].real * b->e[j][0].imag;
      for (k = 1; k < NCOLF; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[j][k].real
                         + a->e[i][k].imag * b->e[j][k].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[j][k].real
                         - a->e[i][k].real * b->e[j][k].imag;
      }
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Matrix multiplication with no adjoints
// c <-- c + a * b
// c <-- c - a * b
// c <-- a * b
// a - funmatrix
// b - funamatrix

void fun_mult_nn_sum(funmatrix *a, funamatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}

void fun_mult_nn_dif(funmatrix *a, funamatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real -= a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag -= a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}

void fun_mult_nn(funmatrix *a, funamatrix *b, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = a->e[i][0].real * b->e[0][j].real
                      - a->e[i][0].imag * b->e[0][j].imag;
      c->e[i][j].imag = a->e[i][0].imag * b->e[0][j].real
                      + a->e[i][0].real * b->e[0][j].imag;
      for (k = 1; k < NCOLF; k++) {
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                         - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                         + a->e[i][k].real * b->e[k][j].imag;
      }
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Copy fun(a)matrices

void fun_mat_copy(funmatrix *src, funmatrix *dest) {
  *dest = *src;
}

void funa_mat_copy(funamatrix *src, funamatrix *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Add result of scalar multiplication on adjoint funmatrix
// c <-- c + s * bdag

void fun_scalar_mult_sum_adj_matrix(funmatrix *b, Real s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += s * b->e[j][i].real;
      c->e[i][j].imag -= s * b->e[j][i].imag;
    }
  }
}

// Add result of scalar multiplication on adjoint funamatrix

void funa_scalar_mult_sum_adj_matrix(funamatrix *b, Real s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real += s * b->e[j][i].real;
      c->e[i][j].imag -= s * b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Add result of scalar multiplication on funmatrix
// c <-- a + s * b
// c <-- c + s * b

void fun_scalar_mult_sum_matrix(funmatrix *b, Real s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void fun_scalar_mult_add_matrix(funmatrix *a, funmatrix *b, Real s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}

// Add result of scalar multiplication on funamatrix

void funa_scalar_mult_sum_matrix(funamatrix *b, Real s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void funa_scalar_mult_add_matrix(funamatrix *a, funamatrix *b, Real s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Scalar multiplication on adjoint of funmatrix
// b <-- s * adag

void fun_scalar_mult_adj_matrix(funmatrix *a, Real s, funamatrix *b) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = s * a->e[j][i].real;
      b->e[i][j].imag = -1.0 * s * a->e[j][i].imag;
    }
  }
}

// Scalar multiplication on adjoint of funamatrix

void funa_scalar_mult_adj_matrix(funamatrix *a, Real s, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      b->e[i][j].real = s * a->e[j][i].real;
      b->e[i][j].imag = -1.0 * s * a->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Scalar multiplication on funmatrix
// b <-- s * a

void fun_scalar_mult_matrix(funmatrix *a, Real s, funmatrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      b->e[i][j].real = s * a->e[i][j].real;
      b->e[i][j].imag = s * a->e[i][j].imag;
    }
  }
}

// Scalar multiplication on funamatrix

void funa_scalar_mult_matrix(funamatrix *a, Real s, funamatrix *b) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = s * a->e[i][j].real;
      b->e[i][j].imag = s * a->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Scaled matrix multiplication with adjoint of first matrix
// c <-- s * c + adag * b
// c <-- s * c - adag * b
// c <-- s * adag * b
// a,b - funamatrices

void funa_scalar_mult_an_sum(funamatrix *a, funamatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real += s * (a->e[k][i].real * b->e[k][j].real
                              + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[k][i].real * b->e[k][j].imag
                              - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}

void funa_scalar_mult_an_dif(funamatrix *a, funamatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real -= s * (a->e[k][i].real * b->e[k][j].real
                              + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag -= s * (a->e[k][i].real * b->e[k][j].imag
                              - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}

void funa_scalar_mult_an(funamatrix *a, funamatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[0][i].real * b->e[0][j].real
                           + a->e[0][i].imag * b->e[0][j].imag);
      c->e[i][j].imag = s * (a->e[0][i].real * b->e[0][j].imag
                           - a->e[0][i].imag * b->e[0][j].real);
      for (k = 1; k < NCOLF; k++) {
        c->e[i][j].real += s * (a->e[k][i].real * b->e[k][j].real
                              + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[k][i].real * b->e[k][j].imag
                              - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Scaled matrix multiplication with adjoint of second matrix
// c <-- c + a * bdag
// c <-- c - a * bdag
// c <-- a * bdag
// a,b - funmatrix

void fun_scalar_mult_na_sum(funmatrix *a, funmatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}

void fun_scalar_mult_na_dif(funmatrix *a, funmatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real -= s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag -= s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}

void fun_scalar_mult_na(funmatrix *a, funmatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[i][0].real * b->e[j][0].real
                           + a->e[i][0].imag * b->e[j][0].imag);
      c->e[i][j].imag = s * (a->e[i][0].imag * b->e[j][0].real
                           - a->e[i][0].real * b->e[j][0].imag);
      for (k = 1; k < NCOLF; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Scaled matrix multiplication with no adjoints
// c <-- c + s * a * b
// c <-- c - s * a * b
// c <-- s * a * b
// a - funmatrix
// b - funamatrix

void fun_scalar_mult_nn_sum(funmatrix *a, funamatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag);
      }
    }
  }
}

void fun_scalar_mult_nn_dif(funmatrix *a, funamatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      for (k = 0; k < NCOLF; k++) {
        c->e[i][j].real -= s * (a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag -= s * (a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag);
      }
    }
  }
}

void fun_scalar_mult_nn(funmatrix *a, funamatrix *b, Real s, matrix *c) {
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[i][0].real * b->e[0][j].real
                           - a->e[i][0].imag * b->e[0][j].imag);
      c->e[i][j].imag = s * (a->e[i][0].imag * b->e[0][j].real
                           + a->e[i][0].real * b->e[0][j].imag);
      for (k = 1; k < NCOLF; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag);
      }
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Subtract result of scalar multiplication on adjoint funmatrix
// c <-- c - s * bdag

void fun_scalar_mult_dif_adj_matrix(funmatrix *b, Real s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= s * b->e[j][i].real;
      c->e[i][j].imag += s * b->e[j][i].imag;
    }
  }
}

// Subtract result of scalar multiplication on adjoint funamatrix

void funa_scalar_mult_dif_adj_matrix(funamatrix *b, Real s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real -= s * b->e[j][i].real;
      c->e[i][j].imag += s * b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Subtract result of scalar multiplication on funmatrix
// c <-- c - s * b

void fun_scalar_mult_dif_matrix(funmatrix *b, Real s, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real -= s * b->e[i][j].real;
      c->e[i][j].imag -= s * b->e[i][j].imag;
    }
  }
}

// Subtract result of scalar multiplication on funamatrix

void funa_scalar_mult_dif_matrix(funamatrix *b, Real s, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= s * b->e[i][j].real;
      c->e[i][j].imag -= s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Subtract two funmatrices
// c <-- c - b
// c <-- a - b

void fun_dif_matrix(funmatrix *b, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real -= b->e[i][j].real;
      c->e[i][j].imag -= b->e[i][j].imag;
    }
  }
}

void fun_sub_matrix(funmatrix *a, funmatrix *b, funmatrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - b->e[i][j].imag;
    }
  }
}

// Subtract two funamatrices

void funa_dif_matrix(funamatrix *b, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[i][j].real;
      c->e[i][j].imag -= b->e[i][j].imag;
    }
  }
}

void funa_sub_matrix(funamatrix *a, funamatrix *b, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - b->e[i][j].imag;
    }
  }
}

void funa_dif_an_matrix(funmatrix *b, funamatrix *c) {
  register int i, j;
  for (i = 0; i < NCOLF; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[j][i].real;
      c->e[i][j].imag += b->e[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Multiply a matrix and a funmatrix 
// c <-- b*c
// c <-- s*b*c
// c <-- b*a

void fun_mat_mult( matrix *b, funmatrix *c) {
  register int i, j, k;
  for(i = 0; i < NCOL; i++)
    for(j = 0; j < NCOLF; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real = b->e[i][k].real * c->e[k][j].real
                              - b->e[i][k].imag * c->e[k][j].imag;
        c->e[i][j].imag = b->e[i][k].imag * c->e[k][j].real
                                + b->e[i][k].real * c->e[k][j].imag;
      }
}

void fun_mat_scalar_mult( matrix *b, Real s, funmatrix *c) {
  register int i, j, k;
  for(i = 0; i < NCOL; i++)
    for(j = 0; j < NCOLF; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real = s * (b->e[i][k].real * c->e[k][j].real
                              - b->e[i][k].imag * c->e[k][j].imag);
        c->e[i][j].imag = s * (b->e[i][k].imag * c->e[k][j].real
                                + b->e[i][k].real * c->e[k][j].imag);
      }
}

void fun_mat_prod(funmatrix *a, matrix *b, funmatrix *c) {
  register int i, j, k;
  for(i = 0; i < NCOL; i++)
    for(j = 0; j < NCOLF; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real = b->e[i][k].real * a->e[k][j].real
                              - b->e[i][k].imag * a->e[k][j].imag;
        c->e[i][j].imag = b->e[i][k].imag * a->e[k][j].real
                                + b->e[i][k].real * a->e[k][j].imag;
      }
}

void fun_mat_prod_dif(funmatrix *a, matrix *b, funmatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOL; i++)
    for(j = 0; j < NCOLF; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real -= b->e[i][k].real * a->e[k][j].real
                              - b->e[i][k].imag * a->e[k][j].imag;
        c->e[i][j].imag -= b->e[i][k].imag * a->e[k][j].real
                                + b->e[i][k].real * a->e[k][j].imag;
      }
}
// -----------------------------------------------------------------

// Multiply a funamatrix and a matrix 
// c <-- c*b
// c <-- s*c*b
// c <-- a*b
// c <-- c + a*b
// c <-- c - a*b

void funa_mat_mult( matrix *b, funamatrix *c) {
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real = c->e[i][k].real * b->e[k][j].real
                              - c->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag = c->e[i][k].imag * b->e[k][j].real
                                + c->e[i][k].real * b->e[k][j].imag;
      }
}

void funa_mat_scalar_mult( matrix *b, Real s, funmatrix *c) {
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real = s * (c->e[i][k].real * b->e[k][j].real
                              - c->e[i][k].imag * b->e[k][j].imag);
        c->e[i][j].imag = s * (c->e[i][k].imag * b->e[k][j].real
                                + c->e[i][k].real * b->e[k][j].imag);
      }
}

void funa_mat_prod( funamatrix *a, matrix *b, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real = a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag = a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag;
      }
}

void funa_mat_prod_sum( funamatrix *a, matrix *b, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real += a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag;
      }
}

void funa_mat_prod_dif( funamatrix *a, matrix *b, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real -= a->e[i][k].real * b->e[k][j].real
                              - a->e[i][k].imag * b->e[k][j].imag;
        c->e[i][j].imag -= a->e[i][k].imag * b->e[k][j].real
                                + a->e[i][k].real * b->e[k][j].imag;
      }
}

void funa_mat_mult_an(funmatrix *a, matrix *b, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++) {
    for(j = 0; j < NCOL; j++) {
      // Init
      c->e[i][j].real = a->e[0][i].real * b->e[0][j].real
                      + a->e[0][i].imag * b->e[0][j].imag;
      c->e[i][j].imag = a->e[0][i].real * b->e[0][j].imag
                      - a->e[0][i].imag * b->e[0][j].real;     
      for( k = 1; k < NCOL; k++){
        c->e[i][j].real = a->e[k][i].real * b->e[k][j].real
                        + a->e[k][i].imag * b->e[k][j].imag;
        c->e[i][j].imag = a->e[k][i].real * b->e[k][j].imag
                        - a->e[k][i].imag * b->e[k][j].real;
      }
    }
  }
}

void funa_scalar_mat_mult_an(funmatrix *a, matrix *b, Real s, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++) {
    for(j = 0; j < NCOL; j++) {
      // Init
      c->e[i][j].real = s * (a->e[0][i].real * b->e[0][j].real
                          + a->e[0][i].imag * b->e[0][j].imag);
      c->e[i][j].imag = s * (a->e[0][i].real * b->e[0][j].imag
                          - a->e[0][i].imag * b->e[0][j].real);     
      for( k = 1; k < NCOL; k++){
        c->e[i][j].real = s * (a->e[k][i].real * b->e[k][j].real
                            + a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag = s * (a->e[k][i].real * b->e[k][j].imag
                            - a->e[k][i].imag * b->e[k][j].real);
      }
    }
  }
}

void funa_mat_prod_an_sum( funmatrix *a, matrix *b, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real += a->e[k][i].real * b->e[k][j].real
                         + a->e[k][i].imag * b->e[k][j].imag;
        c->e[i][j].imag -= a->e[k][i].imag * b->e[k][j].real
                         - a->e[k][i].real * b->e[k][j].imag;
      }
}

void funa_mat_prod_an_dif( funmatrix *a, matrix *b, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real -= a->e[k][i].real * b->e[k][j].real
                         + a->e[k][i].imag * b->e[k][j].imag;
        c->e[i][j].imag += a->e[k][i].imag * b->e[k][j].real
                         - a->e[k][i].real * b->e[k][j].imag;
      }
}

void funa_scalar_mat_prod_an_dif( funmatrix *a, matrix *b, Real s, funamatrix *c){
  register int i, j, k;
  for(i = 0; i < NCOLF; i++)
    for(j = 0; j < NCOL; j++)
      for( k = 0; k < NCOL; k++){
        c->e[i][j].real -= s * (a->e[k][i].real * b->e[k][j].real
                             +  a->e[k][i].imag * b->e[k][j].imag);
        c->e[i][j].imag += s * (a->e[k][i].imag * b->e[k][j].real
                             -  a->e[k][i].real * b->e[k][j].imag);
      }
}

// -----------------------------------------------------------------
// Real trace of funmatrix * adj funmatrix
Real fun_realtrace(funmatrix *a, funmatrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOLF; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Real trace of funamatrix * adj funamatrix
Real funa_realtrace(funamatrix *a, funamatrix *b) {
  register int i, j;
  register Real sum = 0.0;

  for (i = 0; i < NCOLF; i++) {
    for(j = 0; j < NCOL; j++) {
      sum += a->e[i][j].real * b->e[i][j].real
           + a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Complex trace of adj funmatrix adj funmatrix
complex fun_complextrace_an(funmatrix *a, funmatrix *b) {
  register int i, j;
  register Real sumr = 0.0, sumi = 0.0;
  complex sum;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOLF; j++) {
      sumr += a->e[i][j].real * b->e[i][j].real
            + a->e[i][j].imag * b->e[i][j].imag;
      sumi += a->e[i][j].real * b->e[i][j].imag
            - a->e[i][j].imag * b->e[i][j].real;
    }
  }
  sum.real = sumr;
  sum.imag = sumi;
  return sum;
}

// Complex trace of funmatrix * adj funmatrix
complex fun_complextrace_na(funmatrix *a, funmatrix *b) {
  register int i, j;
  register Real sumr = 0.0, sumi = 0.0;
  complex sum;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOLF; j++) {
      sumr += a->e[i][j].real * b->e[i][j].real
            + a->e[i][j].imag * b->e[i][j].imag;
      sumi += a->e[i][j].imag * b->e[i][j].real
            - a->e[i][j].real * b->e[i][j].imag;
    }
  }
  sum.real = sumr;
  sum.imag = sumi;
  return sum;
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Complex trace of adj funamatrix * adj funamatrix
complex funa_complextrace_an(funamatrix *a, funamatrix *b) {
  register int i, j;
  register Real sumr = 0.0, sumi = 0.0;
  complex sum;

  for (i = 0; i < NCOLF; i++) {
    for(j = 0; j < NCOL; j++) {
      sumr += a->e[i][j].real * b->e[i][j].real
            + a->e[i][j].imag * b->e[i][j].imag;
      sumi += a->e[i][j].real * b->e[i][j].imag
            - a->e[i][j].imag * b->e[i][j].real;
    }
  }
  sum.real = sumr;
  sum.imag = sumi;
  return sum;
}

// Complex trace of funmatrix * adj funmatrix
complex funa_complextrace_na(funamatrix *a, funamatrix *b) {
  register int i, j;
  register Real sumr = 0.0, sumi = 0.0;
  complex sum;

  for (i = 0; i < NCOLF; i++) {
    for(j = 0; j < NCOL; j++) {
      sumr += a->e[i][j].real * b->e[i][j].real
            + a->e[i][j].imag * b->e[i][j].imag;
      sumi += a->e[i][j].imag * b->e[i][j].real
            - a->e[i][j].real * b->e[i][j].imag;
    }
  }
  sum.real = sumr;
  sum.imag = sumi;
  return sum;
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Debug functions
#ifdef DEBUGFUNMATRIX
void init_debug_funmatrix(funmatrix *debugFun, funamatrix *debugFuna,
                          matrix *debugMat){
  int i,j;
  for(i=1;i<3;i++){
    for(j=1;j<4;j++){
      debugFuna->e[j-1][i-1].real=i*j;
      debugFuna->e[j-1][i-1].imag=1;
      debugFun->e[i-1][j-1].real=i*j; 
      debugFun->e[i-1][j-1].imag=1; 
    }
    for(j=1;j<3;j++) {debugMat->e[i-1][j-1].real=i*j; debugMat->e[i-1][j-1].imag=1;}
  }

  fun_dumpmat(debugFun);
  funa_dumpmat(debugFuna);
  dumpmat(debugMat);
}

void debug_fun_copy(funmatrix *a, funmatrix *b)
{
  int i,j;
  for(i=0;i<NCOL;i++)
    for(j=0;j<NCOLF;j++){
      b->e[i][j].real = a->e[i][j].real;
      b->e[i][j].imag = a->e[i][j].imag;
    }
}

void debug_funa_copy(funamatrix *a, funamatrix *b)
{
  int i,j;
  for(i=0;i<NCOLF;i++)
    for(j=0;j<NCOL;j++){
      b->e[i][j].real = a->e[i][j].real;
      b->e[i][j].imag = a->e[i][j].imag;
    }
}

void debug_mat_copy(matrix *a, matrix *b)
{
  int i,j;
  for(i=0;i<NCOL;i++)
    for(j=0;j<NCOL;j++){
      b->e[i][j].real = a->e[i][j].real;
      b->e[i][j].imag = a->e[i][j].imag;
    }
}

void debug_funmatrix(){
  funmatrix debugFun;
  funamatrix debugFuna;
  matrix debugMat;
  init_debug_funmatrix(&debugFun,&debugFuna,&debugMat);

  funmatrix f1,g1;
  funamatrix f2,g2;
  matrix f3,g3; //The g's are constants
  complex a;
  a.real = 1;
  a.imag = 1;

  printf("copy and dumpmat test\n");

  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  debug_fun_copy(&debugFun,&g1);
  debug_funa_copy(&debugFuna,&g2);
  debug_mat_copy(&debugMat,&g3);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  dumpmat(&f3);

  // -------
  printf("sum test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_sum_matrix(&g1, &f1);
  funa_sum_matrix(&g2, &f2);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("add test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_add_matrix(&debugFun ,&g1, &f1);
  funa_add_matrix(&debugFuna ,&g2, &f2);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("adjoint test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_adjoint(&g1, &f2);
  funa_adjoint(&g2, &f1);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("clear_mat test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_clear_mat(&f1);
  funa_clear_mat(&f2);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("c_scalar_mult_sum_adj_mat test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_c_scalar_mult_sum_adj_mat(&g1, &a, &f2);
  funa_c_scalar_mult_sum_adj_mat(&g2, &a, &f1);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("c_scalar_mult_dif_adj_mat test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_c_scalar_mult_dif_adj_mat(&g1, &a, &f2);
  funa_c_scalar_mult_dif_adj_mat(&g2, &a, &f1);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("c_scalar_mult_sum_mat test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_c_scalar_mult_sum_mat(&g1, &a, &f1);
  funa_c_scalar_mult_sum_mat(&g2, &a, &f2);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("c_scalar_mult_sum_mat_adj test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_c_scalar_mult_sum_mat_adj(&g1, &a, &f2);
  funa_c_scalar_mult_sum_mat_adj(&g2, &a, &f1);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

// -------
  printf("c_sclar_mult_mat test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_c_scalar_mult_mat(&g1, &a, &f1);
  funa_c_scalar_mult_mat(&g2, &a, &f2);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // -------  

  // -------
  printf("c scalar mult dif mat test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_c_scalar_mult_dif_mat(&g1, &a, &f1);
  funa_c_scalar_mult_dif_mat(&g2, &a, &f2);

  fun_dumpmat(&f1);
  funa_dumpmat(&f2);
  // ------- 

  // ------- 
  printf("funa_mult_an_sum test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  funa_mult_an_sum(&f2, &g2, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("funa_mult_an_dif test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  funa_mult_an_dif(&f2, &g2, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("funa_mult_an test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  funa_mult_an(&f2, &g2, &f3);

  dumpmat(&f3);
  // -------   

// -------
  printf("fun_mult_na_sum test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_mult_na_sum(&f1, &g1, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("fun_mult_na_diff test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_mult_na_dif(&f1, &g1, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("fun_mult_na test \n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_mult_na(&f1, &g1, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("fun_mult_nn_sum test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_mult_nn_sum(&g1, &g2, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("fun_mult_nn_diff test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_mult_nn_dif(&g1, &g2, &f3);

  dumpmat(&f3);
  // -------   

  // -------
  printf("fun_mult_nn test\n");
  debug_fun_copy(&debugFun,&f1);
  debug_funa_copy(&debugFuna,&f2);
  debug_mat_copy(&debugMat,&f3);

  fun_mult_nn(&g1, &g2, &f3);

  dumpmat(&f3);
  // -------   
}

#endif




//TODO: might need functions for the multiplications of the
//type (NCOLF*NCOL) x (NCOL*NCOLF) that give NCOLF*NCOLF matrices
//if so have to add a new struct to deal with them.

//TODO: might need a new type to hold adjoint funmatrices
//alternatively same comment whenever adjoints are to be used
//and use good variable names. might not be a good idea in the long run

//TODO: maybe make a trace function for funmatrix multiplication

