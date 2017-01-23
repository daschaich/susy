// -----------------------------------------------------------------
// Add result of scalar multiplication on matrix
// c <-- a + s * b
// c <-- c + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

void scalar_mult_sum_matrix(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_add_matrix(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_mult_add_matrix(matrix *a, Real s1, matrix *b, Real s2,
		matrix *c) {
	register int i, j;
	for (i = 0; i < NCOL; i++) {
		for (j = 0; j < NCOL; j++) {
			c->e[i][j].real = s1 * a->e[i][j].real + s2 * b->e[i][j].real;
			c->e[i][j].imag = s1 * a->e[i][j].imag + s2 * b->e[i][j].imag;
		}
	}
}
// -----------------------------------------------------------------
