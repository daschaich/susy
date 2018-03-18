// -----------------------------------------------------------------
#include "generic_includes.h"

#define TOLERANCE (0.0001)
#define MAXERRCOUNT 100

Real max_deviation;
double av_deviation;

#define fixsu2(matrix) \
 { \
    CONJG((*matrix).e[0][0],(*matrix).e[1][1]); \
    CONJG((*matrix).e[0][1],(*matrix).e[1][0]); \
    CNEGATE((*matrix).e[1][0],(*matrix).e[1][0]); \
 }

#define fixsu3(matrix) \
 { \
    bj0r = (*matrix).e[0][0].real; \
    bj0i = (*matrix).e[0][0].imag; \
    bj1r = (*matrix).e[0][1].real; \
    bj1i = (*matrix).e[0][1].imag; \
    bj2r = (*matrix).e[0][2].real; \
    bj2i = (*matrix).e[0][2].imag; \
    ar = (*matrix).e[1][2].real; \
    ai = (*matrix).e[1][2].imag; \
    tr = bj1r*ar - bj1i*ai; \
    ti = bj1r*ai + bj1i*ar; \
    ar = (*matrix).e[1][1].real; \
    ai = (*matrix).e[1][1].imag; \
    tr = tr - bj2r*ar + bj2i*ai; \
    ti = ti - bj2r*ai - bj2i*ar; \
    (*matrix).e[2][0].real = tr; \
    (*matrix).e[2][0].imag = -ti; \
    ar = (*matrix).e[1][0].real; \
    ai = (*matrix).e[1][0].imag; \
    tr = bj2r*ar - bj2i*ai; \
    ti = bj2r*ai + bj2i*ar; \
    ar = (*matrix).e[1][2].real; \
    ai = (*matrix).e[1][2].imag; \
    tr = tr - bj0r*ar + bj0i*ai; \
    ti = ti - bj0r*ai - bj0i*ar; \
    (*matrix).e[2][1].real = tr; \
    (*matrix).e[2][1].imag = -ti; \
    ar = (*matrix).e[1][1].real; \
    ai = (*matrix).e[1][1].imag; \
    tr = bj0r*ar - bj0i*ai; \
    ti = bj0r*ai + bj0i*ar; \
    ar = (*matrix).e[1][0].real; \
    ai = (*matrix).e[1][0].imag; \
    tr = tr - bj1r*ar + bj1i*ai; \
    ti = ti - bj1r*ai - bj1i*ar; \
    (*matrix).e[2][2].real = tr; \
    (*matrix).e[2][2].imag = -ti; \
 }

#if NCOL > 3
void fixsu4(matrix *c){
  /* calculate row 3 of SU(4) matrix from other 3 rows, assumed to be already
     orthonormal    */

  register Real bj0r, bj0i, bj1r, bj1i, bj2r, bj2i;
  register Real ar, ai, tr, ti;

  struct {complex e[3][3];} b;  /* 3x3 submatrix        */
  complex sum,t;
  int is,i,j,k,kk;
  is=-1;      /* sign for perms of columns      */
  for(i=0;i<4;i++){   /* loop over location of element of row 3 */
    kk=0;     /* to allow skipping column i in c    */
    for(k=0;k<4;k++){   /* loop on columns of c       */
      if(k!=i){
        for(j=0;j<2;j++)  /* load rows 0,1 of submatrix into b[][]  */
          b.e[j][kk]=c->e[j][k];
        kk++;
      }
    }
    fixsu3((&b));   /* generate row 2 of b        */
    sum.real = sum.imag = 0.;
    kk=0;     /* dot product with row 2 of c      */
    for(k=0;k<4;k++){
      if(k!=i){
        CMUL_J(b.e[2][kk],c->e[2][k],t);
        CSUM(sum,t);
        kk++;
      }
    }
    if(is<0){CNEGATE(sum,sum);} /* minus sign for EVEN columns    */
    c->e[3][i]=sum;
    is=-is;
  }
}
#endif

int check_deviation(Real deviation) {
  if(max_deviation<deviation)
    max_deviation=deviation;
  av_deviation += deviation*deviation;

  if (deviation < TOLERANCE)
    return 0;
  else
    return 1;
}

void reunit_report_problem_matrix(matrix *mat, int i, int dir) {
  int ii,jj;
  union {
    Real fval;
    int ival;
  } ifval;

  printf("Unitarity problem on node %d, site %d, dir %d tolerance=%e\n",
      mynode(),i,dir,TOLERANCE);
  printf("SU(N) matrix:\n");
  for(ii=0;ii<NCOL;ii++){
    for(jj=0;jj<NCOL;jj++){
      printf("%f ",(*mat).e[ii][jj].real);
      printf("%f ",(*mat).e[ii][jj].imag);
    }
    printf("\n");
  }
  printf("repeat in hex:\n");
  for(ii=0;ii<NCOL;ii++){
    for(jj=0;jj<NCOL;jj++){
      ifval.fval = (*mat).e[ii][jj].real;
      printf("%08x ", ifval.ival);
      ifval.fval = (*mat).e[ii][jj].imag;
      printf("%08x ", ifval.ival);
    }
    printf("\n");
  }
  printf("  \n \n");
  fflush(stdout);
}

int reunit(matrix *c) {
  register Real ar;
  register Real c0r, c0i, c1r, c1i;
#if (NCOL>2)
  register Real c2r, c2i, ai;
#if (NCOL==3)
  register Real tr, ti, bj0r, bj0i, bj1r, bj1i, bj2r, bj2i;
#endif
#if (NCOL>3)
  register Real c3r, c3i;
#endif
#endif
  Real deviation;
  int errors;

  errors = 0;
  /* first normalize row 0 */
  ar = (*c).e[0][0].real * (*c).e[0][0].real    /* sum of squares of row */
     + (*c).e[0][0].imag * (*c).e[0][0].imag
     + (*c).e[0][1].real * (*c).e[0][1].real
     + (*c).e[0][1].imag * (*c).e[0][1].imag;
#if (NCOL>2)
  ar += (*c).e[0][2].real * (*c).e[0][2].real
      + (*c).e[0][2].imag * (*c).e[0][2].imag;
#if (NCOL>3)
  ar += (*c).e[0][3].real * (*c).e[0][3].real
      + (*c).e[0][3].imag * (*c).e[0][3].imag;
#endif
#endif

  deviation = fabs(ar - 1.0);
  errors += check_deviation(deviation);

  ar = 1.0 / sqrt( (double)ar);         /* used to normalize row */
  (*c).e[0][0].real *= ar;
  (*c).e[0][0].imag *= ar;
  (*c).e[0][1].real *= ar;
  (*c).e[0][1].imag *= ar;
#if (NCOL>2)
  (*c).e[0][2].real *= ar;
  (*c).e[0][2].imag *= ar;
#if (NCOL>3)
  (*c).e[0][3].real *= ar;
  (*c).e[0][3].imag *= ar;
#endif
#endif

#if (NCOL==2) /* conclude with row 1 */
  /* Save for checking */
  c0r = (*c).e[1][0].real;
  c0i = (*c).e[1][0].imag;
  c1r = (*c).e[1][1].real;
  c1i = (*c).e[1][1].imag;

  fixsu2(c); /* reconstruct row 1 */

  /* check deviation */
  ar = (c0r - (*c).e[1][0].real) * (c0r - (*c).e[1][0].real) +
    (c0i - (*c).e[1][0].imag) * (c0i - (*c).e[1][0].imag) +
    (c1r - (*c).e[1][1].real) * (c1r - (*c).e[1][1].real) +
    (c1i - (*c).e[1][1].imag) * (c1i - (*c).e[1][1].imag);
  deviation = ar;
  errors += check_deviation(deviation);

#else /* NCOL>2 from here to end */

  /* make row 1 orthogonal to row 0 */
  ar = (*c).e[0][0].real * (*c).e[1][0].real +     /* real part of 0 dot 1 */
    (*c).e[0][0].imag * (*c).e[1][0].imag +
    (*c).e[0][1].real * (*c).e[1][1].real +
    (*c).e[0][1].imag * (*c).e[1][1].imag +
    (*c).e[0][2].real * (*c).e[1][2].real +
    (*c).e[0][2].imag * (*c).e[1][2].imag
#if (NCOL>3)
    +(*c).e[0][3].real * (*c).e[1][3].real +
    (*c).e[0][3].imag * (*c).e[1][3].imag
#endif
   ;
  ai = (*c).e[0][0].real * (*c).e[1][0].imag -     /* imag part of 0 dot 1 */
    (*c).e[0][0].imag * (*c).e[1][0].real +
    (*c).e[0][1].real * (*c).e[1][1].imag -
    (*c).e[0][1].imag * (*c).e[1][1].real +
    (*c).e[0][2].real * (*c).e[1][2].imag -
    (*c).e[0][2].imag * (*c).e[1][2].real
#if (NCOL>3)
    +(*c).e[0][3].real * (*c).e[1][3].imag -
    (*c).e[0][3].imag * (*c).e[1][3].real
#endif
   ;

  deviation = ar*ar + ai*ai;
  errors += check_deviation(deviation);

  /* row 1 -= a * row 0 */
  (*c).e[1][0].real -= ar*(*c).e[0][0].real - ai*(*c).e[0][0].imag;
  (*c).e[1][0].imag -= ar*(*c).e[0][0].imag + ai*(*c).e[0][0].real;
  (*c).e[1][1].real -= ar*(*c).e[0][1].real - ai*(*c).e[0][1].imag;
  (*c).e[1][1].imag -= ar*(*c).e[0][1].imag + ai*(*c).e[0][1].real;
  (*c).e[1][2].real -= ar*(*c).e[0][2].real - ai*(*c).e[0][2].imag;
  (*c).e[1][2].imag -= ar*(*c).e[0][2].imag + ai*(*c).e[0][2].real;
#if (NCOL>3)
  (*c).e[1][3].real -= ar*(*c).e[0][3].real - ai*(*c).e[0][3].imag;
  (*c).e[1][3].imag -= ar*(*c).e[0][3].imag + ai*(*c).e[0][3].real;
#endif

  /* normalize row 1 */
  ar = (*c).e[1][0].real * (*c).e[1][0].real    /* sum of squares of row */
     + (*c).e[1][0].imag * (*c).e[1][0].imag
     + (*c).e[1][1].real * (*c).e[1][1].real
     + (*c).e[1][1].imag * (*c).e[1][1].imag
     + (*c).e[1][2].real * (*c).e[1][2].real
     + (*c).e[1][2].imag * (*c).e[1][2].imag;
#if (NCOL>3)
  ar += (*c).e[1][3].real * (*c).e[1][3].real
      + (*c).e[1][3].imag * (*c).e[1][3].imag;
#endif

  deviation = fabs(ar - 1.);
  errors += check_deviation(deviation);

  ar = 1.0 / sqrt( (double)ar);         /* used to normalize row */
  (*c).e[1][0].real *= ar;
  (*c).e[1][0].imag *= ar;
  (*c).e[1][1].real *= ar;
  (*c).e[1][1].imag *= ar;
  (*c).e[1][2].real *= ar;
  (*c).e[1][2].imag *= ar;
#if (NCOL>3)
  (*c).e[1][3].real *= ar;
  (*c).e[1][3].imag *= ar;
#endif

#if (NCOL==3) /* conclude with row 2 */
  /* Save for checking */
  c0r = (*c).e[2][0].real;
  c0i = (*c).e[2][0].imag;
  c1r = (*c).e[2][1].real;
  c1i = (*c).e[2][1].imag;
  c2r = (*c).e[2][2].real;
  c2i = (*c).e[2][2].imag;

  fixsu3(c); /* reconstruct row 2 */

  /* check deviation */
  ar = (c0r - (*c).e[2][0].real) * (c0r - (*c).e[2][0].real) +
    (c0i - (*c).e[2][0].imag) * (c0i - (*c).e[2][0].imag) +
    (c1r - (*c).e[2][1].real) * (c1r - (*c).e[2][1].real) +
    (c1i - (*c).e[2][1].imag) * (c1i - (*c).e[2][1].imag) +
    (c2r - (*c).e[2][2].real) * (c2r - (*c).e[2][2].real) +
    (c2i - (*c).e[2][2].imag) * (c2i - (*c).e[2][2].imag);
  deviation = ar;
  errors += check_deviation(deviation);

#else /* NCOL>3 from here to end */

  /* make row 2 orthogonal to row 0 */
  ar = (*c).e[0][0].real * (*c).e[2][0].real +     /* real part of 0 dot 2 */
    (*c).e[0][0].imag * (*c).e[2][0].imag +
    (*c).e[0][1].real * (*c).e[2][1].real +
    (*c).e[0][1].imag * (*c).e[2][1].imag +
    (*c).e[0][2].real * (*c).e[2][2].real +
    (*c).e[0][2].imag * (*c).e[2][2].imag +
    (*c).e[0][3].real * (*c).e[2][3].real +
    (*c).e[0][3].imag * (*c).e[2][3].imag;
  ai = (*c).e[0][0].real * (*c).e[2][0].imag -     /* imag part of 0 dot 2 */
    (*c).e[0][0].imag * (*c).e[2][0].real +
    (*c).e[0][1].real * (*c).e[2][1].imag -
    (*c).e[0][1].imag * (*c).e[2][1].real +
    (*c).e[0][2].real * (*c).e[2][2].imag -
    (*c).e[0][2].imag * (*c).e[2][2].real +
    (*c).e[0][3].real * (*c).e[2][3].imag -
    (*c).e[0][3].imag * (*c).e[2][3].real;

  deviation = ar*ar + ai*ai;
  errors += check_deviation(deviation);

  /* row 2 -= a * row 0 */
  (*c).e[2][0].real -= ar*(*c).e[0][0].real - ai*(*c).e[0][0].imag;
  (*c).e[2][0].imag -= ar*(*c).e[0][0].imag + ai*(*c).e[0][0].real;
  (*c).e[2][1].real -= ar*(*c).e[0][1].real - ai*(*c).e[0][1].imag;
  (*c).e[2][1].imag -= ar*(*c).e[0][1].imag + ai*(*c).e[0][1].real;
  (*c).e[2][2].real -= ar*(*c).e[0][2].real - ai*(*c).e[0][2].imag;
  (*c).e[2][2].imag -= ar*(*c).e[0][2].imag + ai*(*c).e[0][2].real;
  (*c).e[2][3].real -= ar*(*c).e[0][3].real - ai*(*c).e[0][3].imag;
  (*c).e[2][3].imag -= ar*(*c).e[0][3].imag + ai*(*c).e[0][3].real;

  /* make row 2 orthogonal to row 1 */
  ar = (*c).e[1][0].real * (*c).e[2][0].real +     /* real part of 1 dot 2 */
    (*c).e[1][0].imag * (*c).e[2][0].imag +
    (*c).e[1][1].real * (*c).e[2][1].real +
    (*c).e[1][1].imag * (*c).e[2][1].imag +
    (*c).e[1][2].real * (*c).e[2][2].real +
    (*c).e[1][2].imag * (*c).e[2][2].imag +
    (*c).e[1][3].real * (*c).e[2][3].real +
    (*c).e[1][3].imag * (*c).e[2][3].imag;
  ai = (*c).e[1][0].real * (*c).e[2][0].imag -     /* imag part of 1 dot 2 */
    (*c).e[1][0].imag * (*c).e[2][0].real +
    (*c).e[1][1].real * (*c).e[2][1].imag -
    (*c).e[1][1].imag * (*c).e[2][1].real +
    (*c).e[1][2].real * (*c).e[2][2].imag -
    (*c).e[1][2].imag * (*c).e[2][2].real +
    (*c).e[1][3].real * (*c).e[2][3].imag -
    (*c).e[1][3].imag * (*c).e[2][3].real;

  deviation = ar*ar + ai*ai;
  errors += check_deviation(deviation);

  /* row 2 -= a * row 1 */
  (*c).e[2][0].real -= ar*(*c).e[1][0].real - ai*(*c).e[1][0].imag;
  (*c).e[2][0].imag -= ar*(*c).e[1][0].imag + ai*(*c).e[1][0].real;
  (*c).e[2][1].real -= ar*(*c).e[1][1].real - ai*(*c).e[1][1].imag;
  (*c).e[2][1].imag -= ar*(*c).e[1][1].imag + ai*(*c).e[1][1].real;
  (*c).e[2][2].real -= ar*(*c).e[1][2].real - ai*(*c).e[1][2].imag;
  (*c).e[2][2].imag -= ar*(*c).e[1][2].imag + ai*(*c).e[1][2].real;
  (*c).e[2][3].real -= ar*(*c).e[1][3].real - ai*(*c).e[1][3].imag;
  (*c).e[2][3].imag -= ar*(*c).e[1][3].imag + ai*(*c).e[1][3].real;

  /* normalize row 2 */
  ar = (*c).e[2][0].real * (*c).e[2][0].real +    /* sum of squares of row */
    (*c).e[2][0].imag * (*c).e[2][0].imag +
    (*c).e[2][1].real * (*c).e[2][1].real +
    (*c).e[2][1].imag * (*c).e[2][1].imag +
    (*c).e[2][2].real * (*c).e[2][2].real +
    (*c).e[2][2].imag * (*c).e[2][2].imag +
    (*c).e[2][3].real * (*c).e[2][3].real +
    (*c).e[2][3].imag * (*c).e[2][3].imag;

  deviation = fabs(ar - 1.);
  errors += check_deviation(deviation);

  ar = 1.0 / sqrt( (double)ar);         /* used to normalize row */
  (*c).e[2][0].real *= ar;
  (*c).e[2][0].imag *= ar;
  (*c).e[2][1].real *= ar;
  (*c).e[2][1].imag *= ar;
  (*c).e[2][2].real *= ar;
  (*c).e[2][2].imag *= ar;
  (*c).e[2][3].real *= ar;
  (*c).e[2][3].imag *= ar;

  /* NCOL=4: conclude with row 3 */
  /* Save for checking */
  c0r = (*c).e[3][0].real;
  c0i = (*c).e[3][0].imag;
  c1r = (*c).e[3][1].real;
  c1i = (*c).e[3][1].imag;
  c2r = (*c).e[3][2].real;
  c2i = (*c).e[3][2].imag;
  c3r = (*c).e[3][3].real;
  c3i = (*c).e[3][3].imag;

  fixsu4(c); /* reconstruct row 3 */

  /* check deviation */
  ar = (c0r - (*c).e[3][0].real) * (c0r - (*c).e[3][0].real) +
    (c0i - (*c).e[3][0].imag) * (c0i - (*c).e[3][0].imag) +
    (c1r - (*c).e[3][1].real) * (c1r - (*c).e[3][1].real) +
    (c1i - (*c).e[3][1].imag) * (c1i - (*c).e[3][1].imag) +
    (c2r - (*c).e[3][2].real) * (c2r - (*c).e[3][2].real) +
    (c2i - (*c).e[3][2].imag) * (c2i - (*c).e[3][2].imag) +
    (c3r - (*c).e[3][3].real) * (c3r - (*c).e[3][3].real) +
    (c3i - (*c).e[3][3].imag) * (c3i - (*c).e[3][3].imag);
  deviation = ar;
  errors += check_deviation(deviation);
#endif   /* NCOL */
#endif   /* NCOL */

  return errors;
}

void reunitarize() {
  register matrix *mat;
  register int i,dir;
  register site *s;
  int errcount = 0;
  int errors;

  max_deviation = 0.;
  av_deviation = 0.;

  FORALLSITES(i,s){
    mat = (matrix *)&(s->link);
    errors = reunit( mat );
    errcount += errors;
    if(errors)reunit_report_problem_matrix(mat,i,dir);
    if(errcount > MAXERRCOUNT)
    {
      printf("Unitarity error count exceeded.\n");
      terminate(1);
    }
  }

#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.3e, avrg %.3e\n",
      mynode(), max_deviation, av_deviation);
#endif
  if(max_deviation> TOLERANCE)
  {
    printf("reunitarize: Node %d unitarity problem, maximum deviation=%e\n",
        mynode(),max_deviation);
    errcount++;
    if(errcount > MAXERRCOUNT)
    {
      printf("Unitarity error count exceeded.\n");
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------
