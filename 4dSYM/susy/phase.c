// -----------------------------------------------------------------
// Measure the phase of the pfaffian
// Based on the algorithm in the appendix of hep-lat/0305002
// !!! hep-lat/9903014 provides important definitions and caveats
#include "susy_includes.h"
#define PFA_TOL 1e-12
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return a_i * b_i for two complex vectors (no conjugation!)
double_complex dot(complex *a, complex *b) {
  int i, Ndat = 16 * DIMF;
  complex tc;
  double_complex dot = cmplx(0.0, 0.0);

  for (i = 0; i < sites_on_node * Ndat; i++) {
    CMUL(a[i], b[i], tc);
    CSUM(dot, tc);
  }
  // Accumulate inner product across all nodes
  g_dcomplexsum(&dot);
  return dot;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for the fermion_op matvec
// Convert between complex vectors and Twist_Fermions
// Use plaq_src as temporary storage to be gathered
void matvec(complex *in, complex *out) {
  register site *s;
  int i, j, iter;
  msg_tag *mtag0 = NULL, *mtag1 = NULL, *mtag2;
  su3_vector *plaq23, *plaq13, *plaq12;

  // Copy complex vector into Twist_Fermion src
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  // !!! Need to gather & cycle over fields to ensure non-zero Q[i + 1] M Q[i]
  iter = 0;
  FORALLSITES(i, s) {
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(in[iter]), &(src[i].Fsite.c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Flink[4].c[j]));
      iter++;

      set_complex_equal(&(in[iter]), &(src[i].Flink[0].c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Fplaq[3].c[j]));  // 0, 4 --> 3
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Flink[1].c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Fplaq[6].c[j]));  // 1, 4 --> 6
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Flink[2].c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Fplaq[8].c[j]));  // 2, 4 --> 8
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Flink[3].c[j]));
      iter++;
      set_complex_equal(&(in[iter]), &(src[i].Fplaq[9].c[j]));  // 3, 4 --> 9
      iter++;

      set_complex_equal(&(in[iter]), &(src[i].Fplaq[0].c[j]));  // 0, 1 --> 0
      iter++;
      set_complex_equal(&(in[iter]), &(plaq_src[7][i].c[j]));   // 2, 3 --> 7
      iter++;

      set_complex_equal(&(in[iter]), &(src[i].Fplaq[1].c[j]));  // 0, 2 --> 1
      iter++;
      set_complex_equal(&(in[iter]), &(plaq_src[5][i].c[j]));   // 1, 3 --> 5
      iter++;

      set_complex_equal(&(in[iter]), &(src[i].Fplaq[2].c[j]));  // 0, 3 --> 2
      iter++;
      set_complex_equal(&(in[iter]), &(plaq_src[4][i].c[j]));   // 1, 2 --> 4
      iter++;
    }
  }

  // Gather plaq_src[7] (2, 3) from x - 0 - 1 (gather path 23),
  // plaq_src[5] (1, 3) from x - 0 - 2 (gather path 31)
  // and plaq_src[4] (1, 2) from x - 0 - 3 (gather path 35)
  mtag0 = start_gather_field(plaq_src[7], sizeof(su3_vector),
                             23, EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_field(plaq_src[5], sizeof(su3_vector),
                             31, EVENANDODD, gen_pt[1]);
  mtag2 = start_gather_field(plaq_src[4], sizeof(su3_vector),
                             35, EVENANDODD, gen_pt[2]);

  wait_gather(mtag0);
  wait_gather(mtag1);
  wait_gather(mtag2);
  FORALLSITES(i, s) {
    su3vec_copy((su3_vector *)(gen_pt[0][i]), &(src[i].Fplaq[7]));  // 2, 3
    su3vec_copy((su3_vector *)(gen_pt[1][i]), &(src[i].Fplaq[5]));  // 1, 3
    su3vec_copy((su3_vector *)(gen_pt[2][i]), &(src[i].Fplaq[4]));  // 1, 2
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  cleanup_gather(mtag2);

#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the input vector
  int Ndat = 16 * DIMF;
  if (iter != sites_on_node * Ndat) {
    printf("phase: cycled over %d of %d input components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif

  fermion_op(src, res, PLUS);    // D
//  fermion_op(src, res, MINUS);    // Ddag
  Nmatvecs++;

  // Copy the resulting Twist_Fermion res back to complex vector y
  // Gather plaq_src[7] (2, 3) from x + 0 + 1 (gather path 22),
  // plaq_src[5] (1, 3) from x + 0 + 2 (gather path 30)
  // and plaq_src[4] (1, 2) from x + 0 + 3 (gather path 34)
  FORALLSITES(i, s) {
    su3vec_copy(&(res[i].Fplaq[7]), &(plaq_src[7][i]));
    su3vec_copy(&(res[i].Fplaq[5]), &(plaq_src[5][i]));
    su3vec_copy(&(res[i].Fplaq[4]), &(plaq_src[4][i]));
  }
  mtag0 = start_gather_field(plaq_src[7], sizeof(su3_vector),
                             22, EVENANDODD, gen_pt[0]);
  mtag1 = start_gather_field(plaq_src[5], sizeof(su3_vector),
                             30, EVENANDODD, gen_pt[1]);
  mtag2 = start_gather_field(plaq_src[4], sizeof(su3_vector),
                             34, EVENANDODD, gen_pt[2]);

  wait_gather(mtag0);
  wait_gather(mtag1);
  wait_gather(mtag2);
  iter = 0;
  FORALLSITES(i, s) {
    plaq23 = (su3_vector *)(gen_pt[0][i]);
    plaq13 = (su3_vector *)(gen_pt[1][i]);
    plaq12 = (su3_vector *)(gen_pt[2][i]);
    for (j = 0; j < DIMF; j++) {
      set_complex_equal(&(res[i].Fsite.c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Flink[4].c[j]), &(out[iter]));
      iter++;

      set_complex_equal(&(res[i].Flink[0].c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Fplaq[3].c[j]), &(out[iter])); // 0, 4 --> 3
      iter++;
      set_complex_equal(&(res[i].Flink[1].c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Fplaq[6].c[j]), &(out[iter])); // 1, 4 --> 6
      iter++;
      set_complex_equal(&(res[i].Flink[2].c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Fplaq[8].c[j]), &(out[iter])); // 2, 4 --> 8
      iter++;
      set_complex_equal(&(res[i].Flink[3].c[j]), &(out[iter]));
      iter++;
      set_complex_equal(&(res[i].Fplaq[9].c[j]), &(out[iter])); // 3, 4 --> 9
      iter++;

      set_complex_equal(&(res[i].Fplaq[0].c[j]), &(out[iter])); // 0, 1 --> 0
      iter++;
      set_complex_equal(&(plaq23->c[j]), &(out[iter]));
      iter++;

      set_complex_equal(&(res[i].Fplaq[1].c[j]), &(out[iter])); // 0, 2 --> 1
      iter++;
      set_complex_equal(&(plaq13->c[j]), &(out[iter]));
      iter++;

      set_complex_equal(&(res[i].Fplaq[2].c[j]), &(out[iter])); // 0, 3 --> 2
      iter++;
      set_complex_equal(&(plaq12->c[j]), &(out[iter]));
      iter++;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
  cleanup_gather(mtag2);

#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the output vector
  if (iter != sites_on_node * Ndat) {
    printf("phase: cycled over %d of %d output components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void phase() {
  register int i, j, k;
  int Ndat = 16 * DIMF, shift = this_node * sites_on_node * Ndat;
  double phase, log_mag, tr, dtime;
  complex temp, temp2;
  complex *diag = malloc(volume * Ndat * sizeof(*diag));
  complex *MonC = malloc(sites_on_node * Ndat * sizeof(*MonC));
  complex **Q = malloc(volume * Ndat * sizeof(**Q));

  if (Q == NULL) {
    printf("phase: can't malloc Q\n");
    fflush(stdout);
    terminate(1);
  }

  // Allocate and initialize Q, checking whether or not to reload
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  // We keep part of every column on this node
  // The diagonal elements are distributed across different nodes
  // according to shift = this_node * sites_on_node * Ndat defined above
  // Below we will only set diag[i] on the appropriate node
  // and then scatter it by summing over nodes
  if (ckpt_load < 0)
    ckpt_load = 0;    // Cheap trick to allocate minimum necessary columns
  for (i = ckpt_load; i < volume * Ndat; i++) {
    Q[i] = malloc(sites_on_node * Ndat * sizeof(complex));
    diag[i] = cmplx(0.0, 0.0);    // Initialize to zero
  }
  if (Q[volume * Ndat - 1] == NULL) {
    printf("phase: can't malloc Q[i]\n");
    fflush(stdout);
    terminate(1);
  }

  if (ckpt_load > 0) {    // Load from files
    load_diag(diag, ckpt_load);   // Overwrite initial zeroes above
    loadQ(Q, ckpt_load);
  }
  else {                  // Initialize to zero
    for (i = 0; i < volume * Ndat; i++) {
      for (j = 0; j < sites_on_node * Ndat; j++)
        Q[i][j] = cmplx(0.0, 0.0);
    }
    for (i = 0; i < sites_on_node * Ndat; i++)
      Q[shift + i][i] = cmplx(1.0, 0.0);
  }

  // Cycle over ALL pairs of columns after ckpt_load
  for (i = ckpt_load; i < volume * Ndat - 1; i += 2) {
    // Checkpoint if requested -- make sure ckpt_save = 0 doesn't trigger this
    if (i == ckpt_save || i + 1 == ckpt_save) {
      if (ckpt_save > 0)
        break;
    }

    node0_printf("Columns %d--%d of %d: ", i + 1, i + 2, volume * Ndat);
    Nmatvecs = 0;
    dtime = -dclock();

    // q_{i + 1} --> q_{i + 1} / <q_{i + 1} | D | q_i>
    // But only if <q_{i + 1} | D | q_i> is non-zero!
    matvec(Q[i], MonC);
    temp = dot(Q[i + 1], MonC);
    // !!! Must have non-vanishing matrix element!
    if (cabs_sq(&temp) < PFA_TOL) {
      node0_printf("phase: <%d | D | %d> = (%.4g, %.4g) too small\n",
                   i + 2, i + 1, temp.real, temp.imag);
      terminate(1);
    }

    // All loops over j (the elements of the given column) run in parallel
    // We can shorten many of them since Q is upper-triangular
    if (i + 1 < sites_on_node * Ndat) {
      for (j = 0; j < i + 2; j++) {
        CDIV(Q[i + 1][j], temp, temp2);
        set_complex_equal(&temp2, &(Q[i + 1][j]));
      }
    }
    else {
      for (j = 0; j < sites_on_node * Ndat; j++) {
        CDIV(Q[i + 1][j], temp, temp2);
        set_complex_equal(&temp2, &(Q[i + 1][j]));
      }
    }

    // Cycle over ALL subsequent columns
    for (k = i + 2; k < volume * Ndat; k++) {
      // q_k --> q_k - q_i <q_{i + 1} | D | q_k> - q_{i + 1} <q_i | D | q_k>
      matvec(Q[k], MonC);
      temp = dot(Q[i + 1], MonC);
      if (i + 1 < sites_on_node * Ndat) {
        for (j = 0; j < i + 2; j++) {
          CMUL(Q[i][j], temp, temp2);
          CDIF(Q[k][j], temp2);
        }
      }
      else {
        for (j = 0; j < sites_on_node * Ndat; j++) {
          CMUL(Q[i][j], temp, temp2);
          CDIF(Q[k][j], temp2);
        }
      }

      temp = dot(Q[i], MonC);
      if (i + 1 < sites_on_node * Ndat) {
        for (j = 0; j < i + 2; j++) {
          CMUL(Q[i + 1][j], temp, temp2);
          CSUM(Q[k][j], temp2);
        }
      }
      else {
        for (j = 0; j < sites_on_node * Ndat; j++) {
          CMUL(Q[i + 1][j], temp, temp2);
          CSUM(Q[k][j], temp2);
        }
      }
    }
    // Print some timing information
    // and make sure diagonal elements are still sane
    dtime += dclock();
    node0_printf("%d matvecs in %.4g seconds ", Nmatvecs, dtime);

    // Save diagonal element and free columns i and i+1
    // We use shift to single out the appropriate node, then sum
    // Note that Q[i][i - shift] = 1
    if (0 <= i + 1 - shift && i + 1 - shift < sites_on_node * Ndat)
      set_complex_equal(&(Q[i + 1][i + 1 - shift]), &(diag[i + 1]));
    g_complexsum(&(diag[i + 1]));
    node0_printf("%.16g %.16g\n", diag[i + 1].real, diag[i + 1].imag);
    fflush(stdout);

    // Done with these
    free(Q[i]);
    free(Q[i + 1]);
  }
  if (ckpt_save > 0) {    // Write out from ckpt_save to the end and finish
    save_diag(diag, ckpt_save);
    free(diag);

    saveQ(Q, ckpt_save);
    for (i = ckpt_save; i < volume * Ndat; i++)
      free(Q[i]);
    free(Q);
    free(MonC);
    return;   // Not ready for pfaffian yet
  }
  free(MonC);
  free(Q);

  // Q is triangular by construction
  // Compute its determinant from its non-trivial diagonal elements diag[i + 1]
  // expressed as the phase and (the log of) the magnitude
  // Every node has a copy of the diagonal elements, see above
  phase = 0.0;
  log_mag = 0.0;
  if (this_node == 0) {
    for (i = 1; i < volume * Ndat; i += 2) {   // Start at 1, use diag[i]
      phase += carg(&(diag[i]));
      log_mag += log(cabs_sq(&(diag[i]))) / 2.0;
    }

    // pf(M) = (det Q)^{-1}
    // Negate log of magnitude and phase
    // Following the C++ code, keep phase in [0, 2pi)
    tr = fmod(-1.0 * phase, TWOPI);
    if (tr < 0)
      phase = tr + TWOPI;
    else
      phase = tr;
  }
  node0_printf("PFAFF %.8g %.8g %.8g %.8g\n", -1.0 * log_mag, phase,
               fabs(cos(phase)), fabs(sin(phase)));
  fflush(stdout);
  free(diag);
}
#endif
// -----------------------------------------------------------------
