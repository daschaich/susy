// -----------------------------------------------------------------
// A couple of library-like routines that loop over all sites
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a gauge field as an array of NUMLINK matrices
void copy_bosons(int sign) {
  register int i, j;
  register site *s;
  
  if(sign == 1)
  {
    FORALLSITES(i, s) {
      mat_copy(&s->link, &s->old_link);
      for(j=0;j<NSCALAR;j++)
        mat_copy(&s->X[j], &s->old_X[j]);
    }
  }
  else if(sign == -1)
  {
    FORALLSITES(i, s) {
      mat_copy(&s->old_link, &s->link);
      for(j=0;j<NSCALAR;j++)
        mat_copy(&s->old_X[j], &s->X[j]);
    }
  }
  else {
    node0_printf("Error: incorrect sign in copy_boson: %d\n", sign);
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Shift a matrix without parallel transport
// The dir should come from goffset
void shiftmat(matrix *dat, matrix *temp, int dir) {
  register int i;
  register site *s;
  msg_tag *mtag;

  mtag = start_gather_field(dat, sizeof(matrix),
                            dir, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i, s)
    mat_copy((matrix *)gen_pt[0][i], &(temp[i]));
  cleanup_gather(mtag);
  FORALLSITES(i, s)
    mat_copy(&(temp[i]), &(dat[i]));
}
// -----------------------------------------------------------------
