// -----------------------------------------------------------------
#include "susy_includes.h"

#define MaxJacobiIters 50
#define ROTATE(a, i, j, k, l) (g = a[i][j], \
                               h = a[k][l], \
                               a[i][j] = g - s * (h + g * tau), \
                               a[k][l] = h + s * (g - h * tau));
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Numerical Recipes helper functions
void nrerror(char error_text[]) {
  node0_printf("Numerical Recipes run-time error:\n");
  node0_printf("%s\n", error_text);
  node0_printf("...now exiting to system...\n");
  fflush(stdout);
  exit(1);
}

// Allocate vector with index range v[nl, ..., nh]
Real* vector(int nl, int nh) {
  Real *v = malloc((unsigned)(nh - nl + 1) * sizeof(*v));
  if (!v)
    nrerror("allocation failure in vector()");
  return v - nl;
}

// Allocate matrix with indices range m[nrl, ..., nrh][ncl, ..., .nch]
Real** matrix(int nrl, int nrh, int ncl, int nch) {
  int i;
  Real **m = malloc((unsigned)(nrh - nrl + 1) * sizeof(Real*));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = malloc((unsigned)(nch - ncl + 1) * sizeof(Real));
    if (!m[i])
        nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

// Free vector with index range v[nl, ..., nh]
void free_vector(Real *v, int nl, int nh) {
  free((char*)(v + nl));
}

// Free matrix with indices range m[nrl, ..., nrh][ncl, ..., .nch]
void free_matrix(Real **m, int nrl, int nrh, int ncl, int nch) {
  int i;

  for (i = nrh; i >= nrl; i--)
    free((char*)(m[i] + ncl));
  free((char*)(m + nrl));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Numerical Recipes Jacobi eigenvalue algorithm
void jacobi(Real **a, int n, Real d[], Real **v, int *nrot) {
  int i, j, iq, ip;
  Real tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

  // Initialize
  b = vector(1, n);
  z = vector(1, n);
  for (ip = 1; ip <= n; ip++) {
    for (iq = 1; iq <= n; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }
  for (ip = 1; ip<=n;ip++) {
    b[ip] = a[ip][ip];
    d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  *nrot = 0;
  for (i = 1; i <= MaxJacobiIters; i++) {
    sm = 0.0;
    for (ip = 1; ip <= n - 1; ip++) {
      for (iq = ip + 1; iq <= n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {    // Success -- we're done
      free_vector(z, 1, n);
      free_vector(b, 1, n);
      return;
    }

    if (i < 4)
      tresh = 0.2 * sm / (n * n);
    else
      tresh = 0.0;
    for (ip = 1; ip <= n - 1; ip++) {
      for (iq = ip + 1; iq <=n; iq++) {
        g = 100.0 * fabs(a[ip][iq]);
        if (i > 4 && fabs(d[ip]) + g == fabs(d[ip])
                  && fabs(d[iq]) + g == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h = d[iq] - d[ip];
          if (fabs(h) + g == fabs(h))
            t = (a[ip][iq]) / h;
          else {
            theta = 0.5 * h / (a[ip][iq]);
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / sqrt(1 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq] = 0.0;
          for (j = 1; j <= ip - 1; j++)
            ROTATE(a, j, ip, j, iq)
          for (j = ip + 1; j <= iq - 1; j++)
            ROTATE(a, ip, j, j, iq)
          for (j = iq + 1; j <= n; j++)
            ROTATE(a, ip, j, iq, j)
          for (j = 1; j <= n; j++)
            ROTATE(v, j, ip, j, iq)

          ++(*nrot);
        }
      }
    }
    for (ip = 1; ip <= n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  node0_printf("Jacobi exhausted %d iterations\n", MaxJacobiIters);
  nrerror("Too many iterations in Jacobi routine");
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Given matrix aa, calculate the unitary polar decomposition element
// bb = [1 / (aadag.aa)].u
// Use the Jacobi eigenvalue algorithm from Numerical Recipes
// to diagonalize usquared = aadag.a,
// then project out by the matrix 1 / (aadig.a)
// Not returning the real matrix, but could add that if necessary
void polar(su3_matrix_f *aa, su3_matrix_f *bb) {
  int ci, cf, kt, jt = 1, nmat = 2 * NCOL, ntot;
  Real creal, cimag, testtrace;
  Real **a, **v, *w;
  su3_matrix_f usquared, seig, weig, tmat;

//  node0_printf("enter polar\n");
//  dumpmat_f(aa);

  mult_su3_na_f(aa, aa, &usquared);
//  dumpmat_f(&usquared);

  // Load the input su3_matrix_f aa into a real-valued data structure
  a = matrix(1, nmat, 1, nmat);
  v = matrix(1, nmat, 1, nmat);
  w = vector(1, nmat);
  for (ci = 0; ci < NCOL; ci++) {
    kt = 1;
    for (cf = 0; cf < NCOL; cf++) {
      creal = usquared.e[ci][cf].real;
      cimag = usquared.e[ci][cf].imag;
      a[jt][kt] = creal;
      a[jt + 1][kt] = cimag;
      a[jt][kt + 1] = -cimag;
      a[jt + 1][kt + 1] = creal;
      kt += 2;
    }
    jt += 2;
  }

  // Check input to the Jacobi algorithm
//  printf("Input to Jacobi:\n");
//  for (jt = 1; jt <= nmat; jt++) {
//    for (kt = 1; kt <= nmat; kt++)
//      printf(" %.4g", a[jt][kt]);
//    printf("\n");
//  }

  // The Jacobi routine puts the eigenvectors of a into v,
  // and the eigenvalues of a into w
  jacobi(a, nmat, w, v, &ntot);
//  printf("Jacobi needed %d steps\n", ntot);

  // Move the results back into su3_matrix_f structures
  jt = 1;
  for (ci = 0; ci < NCOL; ci++) {
    kt = 1;
    for (cf = 0; cf < NCOL; cf++) {
      seig.e[ci][cf].real = v[jt][kt];
      seig.e[ci][cf].imag = v[jt + 1][kt];
      weig.e[ci][cf] = cmplx(0.0, 0.0);
      kt += 2;
    }
    jt += 2;
  }

  kt = 1;
  for (cf = 0; cf < NCOL; cf++) {
    weig.e[cf][cf].real = 1.0 / sqrt(w[kt]);
    kt += 2;
  }

  // seig are the eigenvectors
  // weig is 1 / sqrt(w) for eigenvalues w
  mult_su3_na_f(&weig, &seig, &usquared);
  mult_su3_nn_f(&seig, &usquared, &tmat);
  mult_su3_nn_f(&tmat, aa, bb);

  // Print to check the output su3_matrix_f bb
//  printf("check matrix hat u\n");
//  dumpmat_f(bb);

  // Check unitarity of bb
  testtrace = realtrace_su3_f(bb, bb);
  if (fabs(testtrace / (Real)NCOL - 1.0) > 1.0e-6)
    node0_printf("Error getting unitary piece: trace = %.4g\n", testtrace);

  free_vector(w, 1, nmat);
  free_matrix(a, 1, nmat, 1, nmat);
  free_matrix(v, 1, nmat, 1, nmat);
}
#undef ROTATE
// -----------------------------------------------------------------
