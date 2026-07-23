// -----------------------------------------------------------------
// Mostly routines on individual Twist_Fermions,
// which could be moved into the libraries
// The last two are exceptions that loop over all sites
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void dump_TF(Twist_Fermion *in) {
  int mu;
  node0_printf("Fsite:   ");
  dumpmat(&(in->Fsite));
  node0_printf("Fsitez:   ");
  dumpmat(&(in->Fsitez));
  node0_printf("Fsitezb:   ");
  dumpmat(&(in->Fsitezb));
  node0_printf("Fsiteeb:   ");
  dumpmat(&(in->Fsiteeb));
  FORALLDIR(mu) {
    node0_printf("Flink %d: ", mu);
    dumpmat(&(in->Flink[mu]));
    node0_printf("Flinkb %d: ", mu);
    dumpmat(&(in->Flinkb[mu]));
    node0_printf("Fthetab %d: ", mu);
    dumpmat(&(in->Fthetab[mu]));
  }
  for (mu = 0; mu < NPLAQ; mu++) {
    node0_printf("Fplaq %d: ", mu);
    dumpmat(&(in->Fplaq[mu]));
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
void clear_TF(Twist_Fermion *in) {
  register int i;
  clear_mat(&(in->Fsite));
  clear_mat(&(in->Fsitez));
  clear_mat(&(in->Fsitezb));
  clear_mat(&(in->Fsiteeb));
  FORALLDIR(i){
    clear_mat(&(in->Flink[i]));
    clear_mat(&(in->Flinkb[i]));
    clear_mat(&(in->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    clear_mat(&(in->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the squared magnitude of a Twist_Fermion, ReTr[adag.a]
Real magsq_TF(Twist_Fermion *in) {
  register int i;
  register Real sum;
  sum = realtrace(&(in->Fsite), &(in->Fsite));
  sum += realtrace(&(in->Fsitezb), &(in->Fsitezb));
  sum += realtrace(&(in->Fsiteeb), &(in->Fsiteeb));
  sum += realtrace(&(in->Fsitez), &(in->Fsitez));
  FORALLDIR(i){
    sum += realtrace(&(in->Flink[i]), &(in->Flink[i]));
    sum += realtrace(&(in->Flinkb[i]), &(in->Flinkb[i]));
    sum += realtrace(&(in->Fthetab[i]), &(in->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    sum += realtrace(&(in->Fplaq[i]), &(in->Fplaq[i]));
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the dot product of two Twist_Fermions, Tr[adag.b]
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b) {
  register int i;
  complex sum, tc;
  sum = complextrace_an(&(a->Fsite), &(b->Fsite));
  tc = complextrace_an(&(a->Fsitezb), &(b->Fsitezb));
  CSUM(sum, tc);
  tc = complextrace_an(&(a->Fsiteeb), &(b->Fsiteeb));
  CSUM(sum, tc);
  tc = complextrace_an(&(a->Fsitez), &(b->Fsitez));
  CSUM(sum, tc);
  FORALLDIR(i) {
    tc = complextrace_an(&(a->Flink[i]), &(b->Flink[i]));
    CSUM(sum, tc);
    tc = complextrace_an(&(a->Flinkb[i]), &(b->Flinkb[i]));
    CSUM(sum, tc);
    tc = complextrace_an(&(a->Fthetab[i]), &(b->Fthetab[i]));
    CSUM(sum, tc);
  }
  for (i = 0; i < NPLAQ; i++) {
    tc = complextrace_an(&(a->Fplaq[i]), &(b->Fplaq[i]));
    CSUM(sum, tc);
  }
  return sum;
}

// c <-- c + ReTr[adag.b]
void TF_rdot_sum(Twist_Fermion *a, Twist_Fermion *b, Real *c) {
  register int i;
  *c += realtrace(&(a->Fsite), &(b->Fsite));
  *c += realtrace(&(a->Fsitezb), &(b->Fsitezb));
  *c += realtrace(&(a->Fsiteeb), &(b->Fsiteeb));
  *c += realtrace(&(a->Fsitez), &(b->Fsitez));
  FORALLDIR(i){
    *c += realtrace(&(a->Flink[i]), &(b->Flink[i]));
    *c += realtrace(&(a->Flinkb[i]), &(b->Flinkb[i]));
    *c += realtrace(&(a->Fthetab[i]), &(b->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    *c += realtrace(&(a->Fplaq[i]), &(b->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// c <-- c + b
void sum_TF(Twist_Fermion *b, Twist_Fermion *c) {
  register int i;
  sum_matrix(&(b->Fsite), &(c->Fsite));
  sum_matrix(&(b->Fsitezb), &(c->Fsitezb));
  sum_matrix(&(b->Fsiteeb), &(c->Fsiteeb));
  sum_matrix(&(b->Fsitez), &(c->Fsitez));
  FORALLDIR(i){
    sum_matrix(&(b->Flink[i]), &(c->Flink[i]));
    sum_matrix(&(b->Flinkb[i]), &(c->Flinkb[i]));
    sum_matrix(&(b->Fthetab[i]), &(c->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    sum_matrix(&(b->Fplaq[i]), &(c->Fplaq[i]));
}

// c <-- s * b
void scalar_mult_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_matrix(&(b->Fsite), s, &(c->Fsite));
  scalar_mult_matrix(&(b->Fsitezb), s, &(c->Fsitezb));
  scalar_mult_matrix(&(b->Fsiteeb), s, &(c->Fsiteeb));
  scalar_mult_matrix(&(b->Fsitez), s, &(c->Fsitez));
  FORALLDIR(i){
    scalar_mult_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
    scalar_mult_matrix(&(b->Flinkb[i]), s, &(c->Flinkb[i]));
    scalar_mult_matrix(&(b->Fthetab[i]), s, &(c->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_matrix(&(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- c + s * b
void scalar_mult_sum_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_sum_matrix(&(b->Fsite), s, &(c->Fsite));
  scalar_mult_sum_matrix(&(b->Fsitezb), s, &(c->Fsitezb));
  scalar_mult_sum_matrix(&(b->Fsiteeb), s, &(c->Fsiteeb));
  scalar_mult_sum_matrix(&(b->Fsitez), s, &(c->Fsitez));
  FORALLDIR(i){
    scalar_mult_sum_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
    scalar_mult_sum_matrix(&(b->Flinkb[i]), s, &(c->Flinkb[i]));
    scalar_mult_sum_matrix(&(b->Fthetab[i]), s, &(c->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_sum_matrix(&(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- a + s * b
void scalar_mult_add_TF(Twist_Fermion *a, Twist_Fermion *b,
                        Real s, Twist_Fermion *c) {

  register int i;
  scalar_mult_add_matrix(&(a->Fsite), &(b->Fsite), s, &(c->Fsite));
  scalar_mult_add_matrix(&(a->Fsitezb), &(b->Fsitezb), s, &(c->Fsitezb));
  scalar_mult_add_matrix(&(a->Fsiteeb), &(b->Fsiteeb), s, &(c->Fsiteeb));
  scalar_mult_add_matrix(&(a->Fsitez), &(b->Fsitez), s, &(c->Fsitez));
  FORALLDIR(i){
    scalar_mult_add_matrix(&(a->Flink[i]), &(b->Flink[i]), s, &(c->Flink[i]));
    scalar_mult_add_matrix(&(a->Flinkb[i]), &(b->Flinkb[i]), s, &(c->Flinkb[i]));
    scalar_mult_add_matrix(&(a->Fthetab[i]), &(b->Fthetab[i]), s, &(c->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_add_matrix(&(a->Fplaq[i]), &(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- c - s * b
void scalar_mult_dif_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_dif_matrix(&(b->Fsite), s, &(c->Fsite));
  scalar_mult_dif_matrix(&(b->Fsitezb), s, &(c->Fsitezb));
  scalar_mult_dif_matrix(&(b->Fsiteeb), s, &(c->Fsiteeb));
  scalar_mult_dif_matrix(&(b->Fsitez), s, &(c->Fsitez));
  FORALLDIR(i){
    scalar_mult_dif_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
    scalar_mult_dif_matrix(&(b->Flinkb[i]), s, &(c->Flinkb[i]));
    scalar_mult_dif_matrix(&(b->Fthetab[i]), s, &(c->Fthetab[i]));
    }
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_dif_matrix(&(b->Fplaq[i]), s, &(c->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a gauge field as an array of NUMLINK matrices
void gauge_field_copy(field_offset src, field_offset dest) {
  register int i, dir, src2, dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    FORALLDIR(dir) {
      mat_copy((matrix *)F_PT(s, src2), (matrix *)F_PT(s, dest2));
      src2 += sizeof(matrix);
      dest2 += sizeof(matrix);
    }
  }
}

// -- edited 

// Copy a scalar fields
void scalar_field_copy(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  FORALLSITES(i, s) {
    mat_copy((matrix *)F_PT(s, src), (matrix *)F_PT(s, dest));
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
