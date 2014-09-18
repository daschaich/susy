// -----------------------------------------------------------------
// Mostly routines on individual Twist_Fermions,
// which could be moved into the libraries
// The last two are exceptions: they copy or shift for all sites
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void conjTF(Twist_Fermion *src, Twist_Fermion *dest) {
  int mu, nu, j;
  for (j = 0; j < DIMF; j++)
    CONJG(src->Fsite.c[j], dest->Fsite.c[j]);

  for (mu = 0; mu < NUMLINK; mu++) {
    for (j = 0; j < DIMF; j++)
      CONJG(src->Flink[mu].c[j], dest->Flink[mu].c[j]);
  }
  for (mu = 0; mu < NUMLINK; mu++) {
    clearvec(&(dest->Fplaq[mu][mu]));
    for (nu = mu + 1; nu < NUMLINK; nu++) {
      for (j = 0; j < DIMF; j++) {
        CONJG(src->Fplaq[mu][nu].c[j], dest->Fplaq[mu][nu].c[j]);
        CNEGATE(dest->Fplaq[mu][nu].c[j], dest->Fplaq[nu][mu].c[j]);
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void dump_TF(Twist_Fermion *source) {
  int mu, nu;
  node0_printf("Fsite:     ");
  dumpvec(&(source->Fsite));
  for (mu = 0; mu < NUMLINK; mu++) {
    node0_printf("Flink %d:   ", mu);
    dumpvec(&(source->Flink[mu]));
  }
  for (mu = 0; mu < NUMLINK; mu++) {
    for (nu = 0; nu < NUMLINK; nu++) {    // Check for anti-symmetry
//    for (nu = mu + 1; nu < NUMLINK; nu++) {
      node0_printf("Fplaq %d %d: ", mu, nu);
      dumpvec(&(source->Fplaq[mu][nu]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a Twist_Fermion (hardly worth a function)
void copy_TF(Twist_Fermion *src, Twist_Fermion *dest) {
  *dest = *src;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Clear a Twist_Fermion
void clear_TF(Twist_Fermion *vec) {
  register int i, j;
  clearvec(&(vec->Fsite));
  for (i = 0; i < NUMLINK; i++)
    clearvec(&(vec->Flink[i]));
  for (i = 0; i < NUMLINK; i++) {
    for (j = 0; j < NUMLINK; j++)
      clearvec(&(vec->Fplaq[i][j]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the squared magnitude of a Twist_Fermion
Real magsq_TF(Twist_Fermion *vec) {
  register int i, j;
  register Real sum;
  sum = magsq_su3vec(&(vec->Fsite));
  for (i = 0; i < NUMLINK; i++)
    sum += magsq_su3vec(&(vec->Flink[i]));
  for (i = 0; i < NUMLINK; i++) {
    for (j = i + 1; j < NUMLINK; j++)
      sum += magsq_su3vec(&(vec->Fplaq[i][j]));
  }
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the dot product of two Twist_Fermions, adag.b
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b) {
  register int i, j;
  complex temp1, temp2;
  temp1 = su3_dot(&(a->Fsite), &(b->Fsite));
  for (i = 0; i < NUMLINK; i++) {
    temp2 = su3_dot(&(a->Flink[i]), &(b->Flink[i]));
    CSUM(temp1, temp2);
  }
  for (i = 0; i < NUMLINK; i++) {
    for (j = i + 1; j < NUMLINK; j++) {
      temp2 = su3_dot(&(a->Fplaq[i][j]), &(b->Fplaq[i][j]));
      CSUM(temp1, temp2);
    }
  }
  return temp1;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest <-- src1 + s * src2
void scalar_mult_add_TF(Twist_Fermion *src1, Twist_Fermion *src2,
                        Real s, Twist_Fermion *dest) {

  register int i, j, k;
  scalar_mult_add_su3_vector(&(src1->Fsite), &(src2->Fsite),
                             s, &(dest->Fsite));

  for (i = 0; i < NUMLINK; i++) {
    scalar_mult_add_su3_vector(&(src1->Flink[i]), &(src2->Flink[i]),
                               s, &(dest->Flink[i]));
  }
  for (i = 0; i < NUMLINK; i++) {
    clearvec(&(dest->Fplaq[i][i]));
    for (j = i + 1; j < NUMLINK; j++) {
      scalar_mult_add_su3_vector(&(src1->Fplaq[i][j]), &(src2->Fplaq[i][j]),
                                 s, &(dest->Fplaq[i][j]));

      for (k = 0; k < DIMF; k++)
        CNEGATE(dest->Fplaq[i][j].c[k], dest->Fplaq[j][i].c[k]);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void scalar_mult_TF(Twist_Fermion *src, Real s, Twist_Fermion *dest) {
  register int i, j, k;
  scalar_mult_su3_vector(&(src->Fsite), s, &(dest->Fsite));
  for (i = 0; i < NUMLINK; i++)
    scalar_mult_su3_vector(&(src->Flink[i]), s, &(dest->Flink[i]));

  for (i = 0; i < NUMLINK; i++) {
    clearvec(&(dest->Fplaq[i][i]));
    for (j = i + 1; j < NUMLINK; j++) {
      scalar_mult_su3_vector(&(src->Fplaq[i][j]), s, &(dest->Fplaq[i][j]));

      for (k = 0; k < DIMF; k++)
        CNEGATE(dest->Fplaq[i][j].c[k], dest->Fplaq[j][i].c[k]);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a gauge field as an array of NUMLINK su3_matrices
void gauge_field_copy_f(field_offset src, field_offset dest) {
  register int i, dir, src2, dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    for (dir = 0; dir < NUMLINK; dir++) {
      su3mat_copy_f((su3_matrix_f *)F_PT(s, src2),
                    (su3_matrix_f *)F_PT(s, dest2));
      src2 += sizeof(su3_matrix_f);
      dest2 += sizeof(su3_matrix_f);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Shift a matrix without parallel transport
// The tag is the desired content of goffset
void shiftmat(field_offset src, field_offset dest, int tag) {
  register int i;
  register site *s;
  msg_tag *mtag;

  mtag = start_gather_site(src, sizeof(su3_matrix_f),
                           tag, EVENANDODD, gen_pt[0]);
  wait_gather(mtag);
  FORALLSITES(i, s) {
    su3mat_copy_f((su3_matrix_f *)gen_pt[0][i],
                  (su3_matrix_f *)F_PT(s, dest));
  }
  cleanup_gather(mtag);
}
// -----------------------------------------------------------------
