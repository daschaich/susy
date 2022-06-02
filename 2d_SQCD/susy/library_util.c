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
  FORALLDIR(mu) {
    node0_printf("Flink %d: ", mu);
    dumpmat(&(in->Flink[mu]));
  }
  node0_printf("Fplaq:   ");
  dumpmat(&(in->Fplaq));
#ifdef FUNSITE
  node0_printf("FFunsite:   ");
  fun_dumpmat(&(in->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(mu) {
    node0_printf("FFunlink %d: ", mu);
    funa_dumpmat(&(in->Funlink[mu]));
  }
#endif
#ifdef FUNPLAQ
  node0_printf("FFunplaq:   ");
  fun_dumpmat(&(in->Funplaq));
#endif
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
  FORALLDIR(i)
    clear_mat(&(in->Flink[i]));
  clear_mat(&(in->Fplaq));
#ifdef FUNSITE
  fun_clear_mat(&(in->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    funa_clear_mat(&(in->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  fun_clear_mat(&(in->Funplaq));
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the squared magnitude of a Twist_Fermion, ReTr[adag.a]
Real magsq_TF(Twist_Fermion *in) {
  register int i;
  register Real sum;
  sum = realtrace(&(in->Fsite), &(in->Fsite));
  FORALLDIR(i)
    sum += realtrace(&(in->Flink[i]), &(in->Flink[i]));
  sum += realtrace(&(in->Fplaq), &(in->Fplaq));
#ifdef FUNSITE
  sum += fun_realtrace(&(in->Funsite), &(in->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    sum += funa_realtrace(&(in->Funlink[i]),
                          &(in->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  sum += fun_realtrace(&(in->Funplaq), &(in->Funplaq));
#endif
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the dot product of two Twist_Fermions, Tr[adag.b]
complex TF_dot(Twist_Fermion *a, Twist_Fermion *b) {
  register int i;
  complex sum, tc;
  sum = complextrace_an(&(a->Fsite), &(b->Fsite));
  FORALLDIR(i) {
    tc = complextrace_an(&(a->Flink[i]), &(b->Flink[i]));
    CSUM(sum, tc);
  }
  tc = complextrace_an(&(a->Fplaq), &(b->Fplaq));
  CSUM(sum, tc);
#ifdef FUNSITE
  tc = fun_complextrace_an(&(a->Funsite), &(b->Funsite));
  CSUM(sum, tc);
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    tc = funa_complextrace_an(&(a->Funlink[i]),
                              &(b->Funlink[i]));
    CSUM(sum, tc);
  }
#endif
#ifdef FUNPLAQ
  tc = fun_complextrace(&(a->Funplaq), &(b->Funplaq));
  CSUM(sum, tc);
#endif
  return sum;
}

// c <-- c + ReTr[adag.b]
void TF_rdot_sum(Twist_Fermion *a, Twist_Fermion *b, Real *c) {
  register int i;
  *c += realtrace(&(a->Fsite), &(b->Fsite));
  FORALLDIR(i)
    *c += realtrace(&(a->Flink[i]), &(b->Flink[i]));
  *c += realtrace(&(a->Fplaq), &(b->Fplaq));
#ifdef FUNSITE
  *c += fun_realtrace(&(a->Funsite), &(b->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    *c += funa_realtrace(&(a->Funlink[i]),
                         &(b->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  *c += fun_realtrace(&(a->Funplaq), &(b->Funplaq));
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// c <-- c + b
void sum_TF(Twist_Fermion *b, Twist_Fermion *c) {
  register int i;
  sum_matrix(&(b->Fsite), &(c->Fsite));
  FORALLDIR(i)
    sum_matrix(&(b->Flink[i]), &(c->Flink[i]));
  sum_matrix(&(b->Fplaq), &(c->Fplaq));
#ifdef FUNSITE
  fun_sum_matrix(&(b->Funsite), &(c->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    funa_sum_matrix(&(b->Funlink[i]), &(c->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  fun_sum_matrix(&(b->Funplaq), &(c->Funplaq));
#endif
}

// c <-- s * b
void scalar_mult_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_matrix(&(b->Fsite), s, &(c->Fsite));
  FORALLDIR(i)
    scalar_mult_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
  scalar_mult_matrix(&(b->Fplaq), s, &(c->Fplaq));
#ifdef FUNSITE
  fun_scalar_mult_matrix(&(b->Funsite), s, &(c->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    funa_scalar_mult_matrix(&(b->Funlink[i]), s, &(c->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  fun_scalar_mult_matrix(&(b->Funplaq), s, &(c->Funplaq));
#endif
}

// c <-- c + s * b
void scalar_mult_sum_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_sum_matrix(&(b->Fsite), s, &(c->Fsite));
  FORALLDIR(i)
    scalar_mult_sum_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
  scalar_mult_sum_matrix(&(b->Fplaq), s, &(c->Fplaq));
#ifdef FUNSITE
  fun_scalar_mult_sum_matrix(&(b->Funsite), s, &(c->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    funa_scalar_mult_sum_matrix(&(b->Funlink[i]), s, &(c->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  fun_scalar_mult_sum_matrix(&(b->Funplaq), s, &(c->Funplaq));
#endif
}

// c <-- a + s * b
void scalar_mult_add_TF(Twist_Fermion *a, Twist_Fermion *b,
                        Real s, Twist_Fermion *c) {

  register int i;
  scalar_mult_add_matrix(&(a->Fsite), &(b->Fsite), s, &(c->Fsite));
  FORALLDIR(i)
    scalar_mult_add_matrix(&(a->Flink[i]), &(b->Flink[i]), s, &(c->Flink[i]));
  scalar_mult_add_matrix(&(a->Fplaq), &(b->Fplaq), s, &(c->Fplaq));
#ifdef FUNSITE
  fun_scalar_mult_add_matrix(&(a->Funsite), &(b->Funsite), s, &(c->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    funa_scalar_mult_add_matrix(&(a->Funsite),
                                &(b->funlink[i]), s, &(c->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  fun_scalar_mult_add_matrix(&(b->Funplaq), &(b->Funplaq), s, &(c->Funplaq));
#endif
}

// c <-- c - s * b
void scalar_mult_dif_TF(Twist_Fermion *b, Real s, Twist_Fermion *c) {
  register int i;
  scalar_mult_dif_matrix(&(b->Fsite), s, &(c->Fsite));
  FORALLDIR(i)
    scalar_mult_dif_matrix(&(b->Flink[i]), s, &(c->Flink[i]));
  scalar_mult_dif_matrix(&(b->Fplaq), s, &(c->Fplaq));
#ifdef FUNSITE
  fun_scalar_mult_dif_matrix(&(b->Funsite), s, &(c->Funsite));
#endif
#ifdef FUNLINK
  FORALLDIR(i) {
    funa_scalar_mult_dif_matrix(&(b->funlink[i]), s, &(c->Funlink[i]));
  }
#endif
#ifdef FUNPLAQ
  fun_scalar_mult_dif_matrix(&(b->Funplaq), s, &(c->Funplaq));
#endif
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

// Copy a scalar fields
void scalar_field_copy(field_offset src, field_offset dest) {
  register int i;
  register site *s;

  FORALLSITES(i, s) {
    fun_mat_copy((funmatrix *)F_PT(s, src), (funmatrix *)F_PT(s, dest));
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
