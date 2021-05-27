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
  dumpmat(&(in->Fvolume));
  FORALLDIR(mu) {
    node0_printf("Flink %d: ", mu);
    dumpmat(&(in->Flink[mu]));
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
  clear_mat(&(in->Fvolume));
  FORALLDIR(i)
    clear_mat(&(in->Flink[i]));
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
  sum += realtrace(&(in->Fvolume), &(in->Fvolume));
  FORALLDIR(i)
    sum += realtrace(&(in->Flink[i]), &(in->Flink[i]));
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
  tc = complextrace_an(&(a->Fvolume), &(b->Fvolume));
  CSUM(sum, tc);
  FORALLDIR(i) {
    tc = complextrace_an(&(a->Flink[i]), &(b->Flink[i]));
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
  *c += realtrace(&(a->Fvolume), &(b->Fvolume));
  FORALLDIR(i)
    *c += realtrace(&(a->Flink[i]), &(b->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    *c += realtrace(&(a->Fplaq[i]), &(b->Fplaq[i]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// c <-- c + b
void sum_TF(Twist_Fermion *b, Twist_Fermion *c) {
  register int i;
  sum_matrix(&(b->Fsite), &(c->Fsite));
  sum_matrix(&(b->Fvolume), &(c->Fvolume));
  FORALLDIR(i)
    sum_matrix(&(b->Flink[i]), &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    sum_matrix(&(b->Fplaq[i]), &(c->Fplaq[i]));
}

// c <-- c - b
void dif_TF(Twist_Fermion *b, Twist_Fermion *c) {
  register int i;
  dif_matrix(&(b->Fvolume), &(c->Fvolume));
  dif_matrix(&(b->Fsite), &(c->Fsite));
  FORALLDIR(i)
    dif_matrix(&(b->Flink[i]), &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    dif_matrix(&(b->Fplaq[i]), &(c->Fplaq[i]));
}

// c <-- a - b
void sub_TF(Twist_Fermion *a, Twist_Fermion *b, Twist_Fermion *c) {
  register int i;
  sub_matrix(&(a->Fsite), &(b->Fsite), &(c->Fsite));
  sub_matrix(&(a->Fvolume), &(b->Fvolume), &(c->Fvolume));
  FORALLDIR(i)
    sub_matrix(&(a->Flink[i]), &(b->Flink[i]), &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    sub_matrix(&(a->Fplaq[i]), &(b->Fplaq[i]), &(c->Fplaq[i]));
}

// c <-- s * b
void scalar_mult_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_matrix(&(b->Fsite), s, &(c->Fsite));
  scalar_mult_matrix(&(b->Fvolume), s, &(c->Fvolume));
  FORALLDIR(i)
    scalar_mult_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_matrix(&(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- c + s * b
void scalar_mult_sum_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_sum_matrix(&(b->Fsite), s, &(c->Fsite));
  scalar_mult_sum_matrix(&(b->Fvolume), s, &(c->Fvolume));
  FORALLDIR(i)
    scalar_mult_sum_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_sum_matrix(&(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- a + s * b
void scalar_mult_add_TF(Twist_Fermion *a, Twist_Fermion *b,
                        Real s, Twist_Fermion *c) {

  register int i;
  scalar_mult_add_matrix(&(a->Fsite), &(b->Fsite), s, &(c->Fsite));
  scalar_mult_add_matrix(&(a->Fvolume), &(b->Fvolume), s, &(c->Fvolume));
  FORALLDIR(i)
    scalar_mult_add_matrix(&(a->Flink[i]), &(b->Flink[i]), s, &(c->Flink[i]));
  for (i = 0; i < NPLAQ; i++)
    scalar_mult_add_matrix(&(a->Fplaq[i]), &(b->Fplaq[i]), s, &(c->Fplaq[i]));
}

// c <-- c - s * b
void scalar_mult_dif_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_dif_matrix(&(b->Fsite), s, &(c->Fsite));
  scalar_mult_dif_matrix(&(b->Fvolume), s, &(c->Fvolume));
  FORALLDIR(i)
    scalar_mult_dif_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
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
