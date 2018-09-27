// -----------------------------------------------------------------
// Measure the phase of the pfaffian
// Based on the algorithm in the appendix of hep-lat/0305002
// !!! hep-lat/9903014 provides important definitions and caveats
// Note that we work with DIMF=NCOL^2-1 rather than NCOL^2
#include "susy_includes.h"
#define PFA_TOL 1e-12
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return a_i * b_i for two complex vectors (no conjugation!)
double_complex inner(complex *a, complex *b) {
  int i, Ndat = NFERMION * DIMF;
  double_complex dot;

  dot.real = a[0].real * b[0].real - a[0].imag * b[0].imag;
  dot.imag = a[0].imag * b[0].real + a[0].real * b[0].imag;
  for (i = 1; i < sites_on_node * Ndat; i++) {
    dot.real += a[i].real * b[i].real - a[i].imag * b[i].imag;
    dot.imag += a[i].imag * b[i].real + a[i].real * b[i].imag;
  }
  // Accumulate inner product across all nodes
  g_dcomplexsum(&dot);
  return dot;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for the fermion_op matvec
// Convert between complex vectors and fermion matrices
void matvec(complex *in, complex *out) {
  register site *s;
  int i, j, k, iter;

  // Copy complex vector into matrix* src
  // with NFERMION * DIMF non-trivial complex components
  // !!! Need to ensure non-zero Q[i + 1] M Q[i]
  // TODO: Can we rearrange this to avoid all the matrix manipulation?
  iter = 0;
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      clear_mat(&(src[k][i]));
      for (j = 0; j < DIMF; j++) {
        c_scalar_mult_sum_mat(&(Lambda[j]), &(in[iter]), &(src[k][i]));
        iter++;
      }
    }
  }
#ifdef DEBUG_CHECK
  // Check that we didn't miss any components of the input vector
  int Ndat = NFERMION * DIMF;
  if (iter != sites_on_node * Ndat) {
    printf("phase: cycled over %d of %d input components\n",
           iter, sites_on_node * Ndat);
    terminate(1);
  }
#endif

  fermion_op(src, res, PLUS);     // D
//  fermion_op(src, res, MINUS);    // Ddag
  Nmatvecs++;

  // Copy the resulting matrix* res back to complex vector y
  iter = 0;
  FORALLSITES(i, s) {
    for (k = 0; k < NFERMION; k++) {
      for (j = 0; j < DIMF; j++) {
        out[iter] = complextrace_nn(&(res[k][i]), &(Lambda[j]));
        iter++;
      }
    }
  }
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
  int Ndat = NFERMION * DIMF, tot_dat = nt * Ndat;
  int shift = this_node * sites_on_node * Ndat;
  double phase, log_mag, tr, dtime;
  complex tc, tc2;
  complex *diag = malloc(sizeof *diag * tot_dat);
  complex *MonC = malloc(sizeof *MonC * sites_on_node * Ndat);
  complex **Q = malloc(sizeof(complex*) * tot_dat);

  if (Q == NULL) {
    printf("phase: can't malloc Q\n");
    fflush(stdout);
    terminate(1);
  }

  // Allocate and initialize Q, checking whether or not to reload
  // Working with DIMF=NCOL^2-1 rather than NCOL^2
  // We keep part of every column on this node
  // The diagonal elements are distributed across different nodes
  // according to shift = this_node * sites_on_node * Ndat defined above
  // Below we will only set diag[i] on the appropriate node
  // and then scatter it by summing over nodes
  if (ckpt_load < 0)
    ckpt_load = 0;    // Cheap trick to allocate minimum necessary columns
  for (i = ckpt_load; i < tot_dat; i++) {
    Q[i] = malloc(sizeof(complex) * sites_on_node * Ndat);
    diag[i] = cmplx(0.0, 0.0);    // Initialize to zero
  }
  if (Q[tot_dat - 1] == NULL) {
    printf("phase: can't malloc Q[i]\n");
    fflush(stdout);
    terminate(1);
  }

  if (ckpt_load > 0) {    // Load from files
    load_diag(diag, ckpt_load);   // Overwrite initial zeroes above
    loadQ(Q, ckpt_load);
  }
  else {                  // Initialize to unit matrix
    for (i = 0; i < tot_dat; i++) {
      for (j = 0; j < sites_on_node * Ndat; j++)
        Q[i][j] = cmplx(0.0, 0.0);
    }
    for (i = 0; i < sites_on_node * Ndat; i++)
      Q[shift + i][i] = cmplx(1.0, 0.0);
  }

  // Cycle over ALL pairs of columns after ckpt_load
  for (i = ckpt_load; i < tot_dat - 1; i += 2) {
    // Checkpoint if requested -- make sure ckpt_save = 0 doesn't trigger this
    if (i == ckpt_save || i + 1 == ckpt_save) {
      if (ckpt_save > 0)
        break;
    }

    node0_printf("Columns %d--%d of %d: ", i + 1, i + 2, tot_dat);
    Nmatvecs = 0;
    dtime = -dclock();

    // q_{i + 1} --> q_{i + 1} / <q_{i + 1} | D | q_i>
    // But only if <q_{i + 1} | D | q_i> is non-zero!
    matvec(Q[i], MonC);
    tc = inner(Q[i + 1], MonC);
    // !!! Must have non-vanishing matrix element!
    if (cabs_sq(&tc) < PFA_TOL) {
      node0_printf("phase: <%d | D | %d> = (%.4g, %.4g) too small\n",
                   i + 2, i + 1, tc.real, tc.imag);
      terminate(1);
    }

    // All loops over j (the elements of the given column) run in parallel
    // We can shorten many of them since Q is upper-triangular
    if (i + 1 < sites_on_node * Ndat) {
      for (j = 0; j < i + 2; j++) {
        CDIV(Q[i + 1][j], tc, tc2);
        set_complex_equal(&tc2, &(Q[i + 1][j]));
      }
    }
    else {
      for (j = 0; j < sites_on_node * Ndat; j++) {
        CDIV(Q[i + 1][j], tc, tc2);
        set_complex_equal(&tc2, &(Q[i + 1][j]));
      }
    }

    // Cycle over ALL subsequent columns
    for (k = i + 2; k < tot_dat; k++) {
      // q_k --> q_k - q_i <q_{i + 1} | D | q_k> - q_{i + 1} <q_i | D | q_k>
      matvec(Q[k], MonC);
      tc = inner(Q[i + 1], MonC);
      if (i + 1 < sites_on_node * Ndat) {
        for (j = 0; j < i + 2; j++)
          CMULDIF(Q[i][j], tc, Q[k][j]);
      }
      else {
        for (j = 0; j < sites_on_node * Ndat; j++)
          CMULDIF(Q[i][j], tc, Q[k][j]);
      }

      tc = inner(Q[i], MonC);
      if (i + 1 < sites_on_node * Ndat) {
        for (j = 0; j < i + 2; j++)
          CMULSUM(Q[i + 1][j], tc, Q[k][j]);
      }
      else {
        for (j = 0; j < sites_on_node * Ndat; j++)
          CMULSUM(Q[i + 1][j], tc, Q[k][j]);
      }
    }
    // Print some timing information
    // and make sure diagonal elements are still sane
    dtime += dclock();
    node0_printf("%li matvecs in %.4g seconds ", Nmatvecs, dtime);

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
    for (i = ckpt_save; i < tot_dat; i++)
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
    for (i = 1; i < tot_dat; i += 2) {   // Start at 1, use diag[i]
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
