// -----------------------------------------------------------------
#include "susy_includes.h"

static Real at, bt, ct, maxarg1, maxarg2;
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(a, b) (maxarg1 = (a), maxarg2 = (b), \
                   maxarg1  >  maxarg2  ?  maxarg1 : maxarg2)
#define PYTHAG(a, b) ((at = fabs(a)) > (bt = fabs(b)) ? \
                      (ct = bt / at, at * sqrt(1.0 + ct * ct)) : \
                      (bt ? (ct = at / bt, bt * sqrt(1.0 + ct * ct)) : 0.0))
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Numerical Recipes singular value decomposition helper functions
void nrerror(char error_text[]) {
  node0_printf("Numerical Recipes run-time error:\n");
  node0_printf("%s\n", error_text);
  node0_printf("...now exiting to system...\n");
  exit(1);
}

Real* vector(int nl, int nh) {
  Real *v = malloc((unsigned)(nh - nl + 1) * sizeof(*v));
  if (!v)
    nrerror("allocation failure in vector()");
  return v - nl;
}

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

void free_vector(Real *v, int nl, int nh) {
  free((char*)(v + nl));
}

void free_matrix(Real **m, int nrl, int nrh, int ncl, int nch) {
  int i;

  for (i = nrh; i >= nrl; i--)
    free((char*)(m[i] + ncl));
  free((char*)(m + nrl));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Numerical Recipes singular value decomposition
void svdcmp(Real **a, int m, int n, Real *w, Real **v) {
  int flag, i, its, j, jj, k, l, nm;
  Real c, f, h, s, x, y, z;
  Real anorm = 0.0, g = 0.0, scale = 0.0;
  Real *rv1, *vector();

  if (m < n)
    nrerror("SVDCMP: You must augment A with extra zero rows");
  rv1 = vector(1, n);
  for (i = 1; i <= n; i++) {
    l = i + 1;
    rv1[i] = scale*g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f = a[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        if (i != n) {
          for (j = l;j<= n;j++) {
            for (s = 0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
            f = s/h;
            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
          }
        }
        for (k = i; k <= m; k++)
          a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i <= m && i != n) {
      for (k = l;k<= n;k++)
        scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<= n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f = a[i][l];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<= n;k++) rv1[k] = a[i][k]/h;
        if (i != m) {
          for (j = l;j<=m;j++) {
            for (s = 0.0,k=l;k<= n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<= n;k++) a[j][k] += s*rv1[k];
          }
        }
        for (k=l;k<= n;k++) a[i][k] *= scale;
      }
    }
    anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }
  for (i = n;i>=1;i--) {
    if (i < n) {
      if (g) {
        for (j = l;j<= n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j = l;j<= n;j++) {
          for (s = 0.0,k=l;k<= n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<= n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j = l;j<= n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g = rv1[i];
    l=i;
  }
  for (i = n;i>=1;i--) {
    l=i+1;
    g = w[i];
    if (i < n)
      for (j = l;j<= n;j++) a[i][j]=0.0;
    if (g) {
      g = 1.0/g;
      if (i != n) {
        for (j = l;j<= n;j++) {
          for (s = 0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
          f = (s/a[i][i])*g;
          for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
        }
      }
      for (j = i;j<=m;j++) a[j][i] *= g;
    } else {
      for (j = i;j<=m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }
  for (k = n;k>=1;k--) {
    for (its = 1;its<=30;its++) {
      flag = 1;
      for (l=k;l>=1;l--) {
        nm=l-1;
        if (fabs(rv1[l])+anorm == anorm) {
          flag = 0;
          break;
        }
        if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l;i<=k;i++) {
          f = s*rv1[i];
          if (fabs(f)+anorm != anorm) {
            g = w[i];
            h=PYTHAG(f,g);
            w[i]=h;
            h=1.0/h;
            c = g*h;
            s = (-f*h);
            for (j = 1;j<=m;j++) {
              y = a[j][nm];
              z = a[j][i];
              a[j][nm] = y*c+z*s;
              a[j][i] = z*c-y*s;
            }
          }
        }
      }
      z = w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j = 1;j<= n;j++) v[j][k]=(-v[j][k]);
        }
        break;
      }
      if (its == 30)
        nrerror("No convergence in 30 SVDCMP iterations");
      x = w[l];
      nm=k-1;
      y = w[nm];
      g = rv1[nm];
      h=rv1[k];
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g = PYTHAG(f,1.0);
      f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c = s = 1.0;
      for (j = l;j<= nm;j++) {
        i = j+1;
        g = rv1[i];
        y = w[i];
        h = s*g;
        g = c*g;
        z = PYTHAG(f,h);
        rv1[j] = z;
        c = f/z;
        s = h/z;
        f = x*c+g*s;
        g = g*c-x*s;
        h = y*s;
        y = y*c;
        for (jj = 1;jj<= n;jj++) {
          x=v[jj][j];
          z = v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i] = z*c-x*s;
        }
        z = PYTHAG(f,h);
        w[j] = z;
        if (z) {
          z = 1.0/z;
          c = f*z;
          s = h*z;
        }
        f = (c*g)+(s*y);
        x=(c*y)-(s*g);
        for (jj = 1;jj<=m;jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y*c+z*s;
          a[jj][i] = z*c-y*s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free_vector(rv1, 1, n);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Given a matrix aa calculate the polar decomposition elements:
// the unitary matrix bb and the remainder cc
// If aa = u w v^T where w is diagonal, then b = u v^dag and b= v w v^dag
// The code is a wrapper around Numerical Recipes routines
void polar(su3_matrix_f *aa, su3_matrix_f *bb, su3_matrix_f *cc) {
  int ci, cf, kt, jt = 1, it, nmat = 2 * NCOL;
  Real creal, cimag;
  Real **a, **v, *w, **p, **u;
  Real testtrace;

//  node0_printf("enter polar\n");
//  dumpmat_f(aa);

  // Load the input su3_matrix_f aa into a real-valued data structure
  a = matrix(1, nmat, 1, nmat);
  v = matrix(1, nmat, 1, nmat);
  w = vector(1, nmat);
  for (ci = 0; ci < NCOL; ci++) {
    kt = 1;
    for (cf = 0; cf < NCOL; cf++) {
      creal = aa->e[ci][cf].real;
      cimag = aa->e[ci][cf].imag;
      a[jt][kt] = creal;
      a[jt + 1][kt] = cimag;
      a[jt][kt + 1] = -cimag;
      a[jt + 1][kt + 1] = creal;
      kt += 2;
    }
    jt += 2;
  }

  // Check input to the SVD
//  printf("Input to SVD:\n");
//  for (jt = 1; jt <= nmat; jt++) {
//    for (kt = 1; kt <= nmat; kt++)
//      printf(" %.4g", a[jt][kt]);
//    printf("\n");
//  }

  // Call the SVD
  // a = U W V^T is overwritten by U
  svdcmp(a, nmat, nmat, w, v);

  // Construct u = a v^T and P = v W v^T
  u = matrix(1, nmat, 1, nmat);
  for (jt = 1; jt <= nmat; jt++) {
    for (kt = 1; kt <= nmat; kt++) {
      u[jt][kt] = a[jt][1] * v[kt][1];          // Initialize
      for (it = 2; it <= nmat; it++)
        u[jt][kt] += a[jt][it] * v[kt][it];
    }
  }

  p = matrix(1, nmat, 1, nmat);
  for (jt = 1; jt <= nmat; jt++) {
    for (kt = 1; kt <= nmat; kt++) {
      p[jt][kt] = v[jt][1] * w[1] * v[kt][1];   // Initialize
      for (it = 2; it <= nmat; it++)
        p[jt][kt] += v[jt][it] * w[it] * v[kt][it];
    }
  }

  // Check the polar decomposition
//  Real **testmat;
//  testmat = matrix(1, nmat, 1, nmat);
//  for (jt =1; jt <= nmat; jt++) {
//    for (kt = 1; kt <= nmat; kt++)
//      testmat[jt][kt] = 0.0;
//  }
//  for (jt = 1; jt <= nmat; jt++) {
//    for (kt = 1; kt <= nmat; kt++) {
//      for (it = 1; it <= nmat; it++)
//        testmat[jt][kt] += u[jt][it] * p[it][kt];
//    }
//  }
//  printf("Test of polar decomposition:\n");
//  for (jt = 1; jt <= nmat; jt++) {
//    for (kt =1; kt <= nmat; kt++)
//      printf(" %.4g", testmat[jt][kt]);
//    printf("\n");
//  }

  // Finally, move the results into the su3_matrix_f output matrices
  jt = 1;
  for (ci = 0; ci < NCOL; ci++) {
    kt = 1;
    for (cf = 0; cf < NCOL; cf++) {
      bb->e[ci][cf].real = u[jt][kt];
      bb->e[ci][cf].imag = u[jt + 1][kt];
      kt += 2;
    }
    jt += 2;
  }

  jt = 1;
  for (ci = 0; ci < NCOL; ci++) {
    kt = 1;
    for (cf = 0; cf < NCOL; cf++) {
      cc->e[ci][cf].real = p[jt][kt];
      cc->e[ci][cf].imag = p[jt + 1][kt];
      kt += 2;
    }
    jt += 2;
  }

  // Print to check that u and p were reloaded properly
//  printf("check reload of matrix u\n");
//  dumpmat_f(bb);
//  printf("check reload of matrix p\n");
//  dumpmat_f(cc);

  // Check unitarity of bb
  testtrace = realtrace_su3_f(bb, bb);
  if (fabs(testtrace / (Real)NCOL - 1.0) > 1.e-6)
    node0_printf("Error getting unitary piece\n");

  free_vector(w, 1, nmat);
  free_matrix(a, 1, nmat, 1, nmat);
  free_matrix(v, 1, nmat, 1, nmat);
  free_matrix(u, 1, nmat, 1, nmat);
  free_matrix(p, 1, nmat, 1, nmat);
}
#undef SIGN
#undef MAX
#undef PYTHAG
// -----------------------------------------------------------------
