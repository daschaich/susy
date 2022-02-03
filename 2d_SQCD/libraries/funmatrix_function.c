// A file that containts the functions needed for funmatrix operations.
// TO BE split in different files at the end.

//Matrix multiplication of funmatrix * (funmatrix)^dag
// c <-- a * b^dag
void fun_scalar_mult(funmatrix *a, funmatrix *b, Real s, matrix *c)
{
  register int i, j, k;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOLF; j++) {
      // Initialize
      c->e[i][j].real = s * (a->e[i][0].real * b->e[j][0].real
                           + a->e[i][0].imag * b->e[j][0].imag);
      c->e[i][j].imag = s * (a->e[i][0].imag * b->e[j][0].real
                           - a->e[i][0].real * b->e[j][0].imag);
      for (k = 1; k < NCOL; k++) {
        c->e[i][j].real += s * (a->e[i][k].real * b->e[j][k].real
                              + a->e[i][k].imag * b->e[j][k].imag);
        c->e[i][j].imag += s * (a->e[i][k].imag * b->e[j][k].real
                              - a->e[i][k].real * b->e[j][k].imag);
      }
    }
  }
}