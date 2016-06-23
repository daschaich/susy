// -----------------------------------------------------------------
// Eigenvalue computation and helper functions
// !!!Path to primme.h may need local customization
#include "susy_includes.h"
#include "../PRIMME/PRIMMESRC/COMMONSRC/primme.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gaussian random Twist_Fermion
void rand_TFsource(Twist_Fermion *src) {
  register int i, j, mu;
  register site *s;
  complex grn;

  // Begin with pure gaussian random numbers
  FORALLSITES(i, s) {
    clear_TF(&(src[i]));
    for (j = 0; j < DIMF; j++) {                // Site fermions
#ifdef SITERAND
      grn.real = gaussian_rand_no(&(s->site_prn));
      grn.imag = gaussian_rand_no(&(s->site_prn));
#else
      grn.real = gaussian_rand_no(&node_prn);
      grn.imag = gaussian_rand_no(&node_prn);
#endif
      c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fsite));
      FORALLDIR(mu) {                           // Link fermions
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Flink[mu]));
      }
      for (mu = 0; mu < NPLAQ; mu++) {         // Plaquette fermions
#ifdef SITERAND
        grn.real = gaussian_rand_no(&(s->site_prn));
        grn.imag = gaussian_rand_no(&(s->site_prn));
#else
        grn.real = gaussian_rand_no(&node_prn);
        grn.imag = gaussian_rand_no(&node_prn);
#endif
        c_scalar_mult_sum_mat(&(Lambda[j]), &grn, &(src[i].Fplaq[mu]));
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get Nvec vectors (stored consecutively) and hit them by the matrix
void av_ov (void *x, void *y, int *Nvec, primme_params *primme) {
  register site *s;
  int i, j, k, mu, iter, ivec, Ndat = 16 * NCOL * NCOL;
  Complex_Z *xx;

  for (ivec = 0; ivec < *Nvec; ivec++) {
    // Copy double precision complex vector x
    // into Real precision Twist_Fermion src
    // Each Twist_Fermion has Ndat=16DIMF non-trivial complex components
    xx = ((Complex_Z*) x) + Ndat * ivec * sites_on_node;  // This vector in x
    iter = 0;
    FORALLSITES(i, s) {
      for (j = 0; j < NCOL; j++) {
        for (k = 0; k < NCOL; k++) {
          src[i].Fsite.e[j][k].real = xx[iter].r;
          src[i].Fsite.e[j][k].imag = xx[iter].i;
          iter++;
          FORALLDIR(mu) {
            src[i].Flink[mu].e[j][k].real = xx[iter].r;
            src[i].Flink[mu].e[j][k].imag = xx[iter].i;
            iter++;
          }
          for (mu = 0; mu < NPLAQ; mu++) {
            src[i].Fplaq[mu].e[j][k].real = xx[iter].r;
            src[i].Fplaq[mu].e[j][k].imag = xx[iter].i;
            iter++;
          }
        }
      }
    }

#ifdef DEBUG_CHECK
    if (iter != Ndat * sites_on_node)
      printf("av_ov: iter = %d after source\n", iter);

    // Check that src is the same magnitude as x[ivec]
    register site *s;
    double xmag = 0.0, src_mag = 0.0;
    xx = ((Complex_Z*) x) + Ndat * ivec * sites_on_node;   // This vector in x
    for (i = 0; i < sites_on_node * Ndat; i++)
      xmag += xx[i].r * xx[i].r + xx[i].i * xx[i].i;
    FORALLSITES(i, s)
      src_mag += magsq_TF(&(src[i]));
    if (fabs(xmag - src_mag) > eig_tol * eig_tol) {
      node0_printf("av_ov: |x[%d]|^2 = %.4g but |src|^2 = %.4g (%.4g)\n",
                   ivec, xmag, src_mag, fabs(xmag - src_mag));
    }
#endif

#ifdef DEBUG_CHECK
    // Check that src is being copied appropriately
    node0_printf("eigVec[0] copy check:\n");
    dump_TF(&(src[0]));
#endif

    DSq(src, res);    // D^2 + fmass^2

    // Copy the resulting Twist_Fermion res back to complex vector y
    // Each Twist_Fermion has Ndat=16DIMF non-trivial complex components
    xx = ((Complex_Z*) y) + Ndat * ivec * sites_on_node;  // This vector in y
    iter = 0;
    FORALLSITES(i, s) {
      for (j = 0; j < NCOL; j++) {
        for (k = 0; k < NCOL; k++) {
          xx[iter].r = (double)res[i].Fsite.e[j][k].real;
          xx[iter].i = (double)res[i].Fsite.e[j][k].imag;
          iter++;
          FORALLDIR(mu) {
            xx[iter].r = (double)res[i].Flink[mu].e[j][k].real;
            xx[iter].i = (double)res[i].Flink[mu].e[j][k].imag;
            iter++;
          }
          for (mu = 0; mu < NPLAQ; mu++) {
            xx[iter].r = (double)res[i].Fplaq[mu].e[j][k].real;
            xx[iter].i = (double)res[i].Fplaq[mu].e[j][k].imag;
            iter++;
          }
        }
      }
    }

#ifdef DEBUG_CHECK
    if (iter != Ndat * sites_on_node)
      printf("av_ov: iter = %d after source\n", iter);

    // Check that res is the same magnitude as y[ivec]
    double ymag = 0.0, res_mag = 0.0;
    xx = ((Complex_Z*) y) + Ndat * ivec * sites_on_node;   // This vector in x
    for (i = 0; i < sites_on_node * Ndat; i++)
      ymag += xx[i].r * xx[i].r + xx[i].i * xx[i].i;
    FORALLSITES(i, s)
      res_mag += magsq_TF(&(res[i]));
    if (fabs(ymag - res_mag) > eig_tol * eig_tol) {
      node0_printf("av_ov: |y[%d]|^2 = %.4g but |res|^2 = %.4g (%.4g)\n",
                   ivec, ymag, res_mag, fabs(ymag - res_mag));
    }
#endif
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function par_GlobalSumDouble is set as primme.globalSumDouble
void par_GlobalSumDouble(void *sendBuf, void *recvBuf,
                         int *count, primme_params *primme) {

  int i;
  for (i = 0; i < *count; i++)
    *((double*)recvBuf + i) = *((double*)sendBuf + i);

  g_vecdoublesum((double*)recvBuf, *count);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Function make_evs computes eigenvalues through PRIMME
// Prints them with a quick check of |D^dag D phi - lambda phi|^2
// If flag==1 we calculate the smallest eigenvalues
// If flag==-1 we calculate the largest eigenvalues
int make_evs(int Nvec, Twist_Fermion **eigVec, double *eigVal, int flag) {
  register site* s;
  int i, j, k, mu, ivec, iter = 0, ret, Ndat = 16 * NCOL * NCOL;
  int maxn = sites_on_node * Ndat;
  double check, *rnorms = malloc(Nvec * sizeof(*rnorms));
  Complex_Z *workVecs = malloc(Nvec * maxn * sizeof(*workVecs));
  static primme_params primme;
  Twist_Fermion tTF, *tmpTF = malloc(sites_on_node * sizeof(*tmpTF));

  // Check memory allocations
  if (workVecs == NULL) {
    node0_printf("ERROR in make_evs: couldn't allocate workVecs\n");
    exit(1);
  }
  if (rnorms == NULL) {
    node0_printf("ERROR in make_evs: couldn't allocate rnorms\n");
    exit(1);
  }

  // Initialize all the eigenvectors to random vectors
  for (ivec = 0; ivec < Nvec; ivec++) {
    eigVal[ivec] = 1e16;
    rand_TFsource(eigVec[ivec]);
  }

  // Copy initial guesses into double-precision temporary fields
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  for (ivec = 0; ivec < Nvec; ivec++) {
    iter = Ndat * ivec * sites_on_node;   // This vector in workvecs
    FORALLSITES(i, s) {
      for (j = 0; j < NCOL; j++) {
        for (k = 0; k < NCOL; k++) {
          workVecs[iter].r = eigVec[ivec][i].Fsite.e[j][k].real;
          workVecs[iter].i = eigVec[ivec][i].Fsite.e[j][k].imag;
          iter++;
          FORALLDIR(mu) {
            workVecs[iter].r = eigVec[ivec][i].Flink[mu].e[j][k].real;
            workVecs[iter].i = eigVec[ivec][i].Flink[mu].e[j][k].imag;
            iter++;
          }
          for (mu = 0; mu < NPLAQ; mu++) {
            workVecs[iter].r = eigVec[ivec][i].Fplaq[mu].e[j][k].real;
            workVecs[iter].i = eigVec[ivec][i].Fplaq[mu].e[j][k].imag;
            iter++;
          }
        }
      }
    }
  }
#ifdef DEBUG_CHECK
  if (iter != Nvec * maxn)
    printf("make_evs: iter = %d after input\n", iter);
#endif

  // Set the parameters of the EV finder
  primme_initialize(&primme);
  primme.n = maxn * number_of_nodes;              // Global size of matrix
  primme.nLocal = maxn;                           // Local volume
  primme.maxOuterIterations = maxIter;
  primme.maxMatvecs = maxIter + 5;
  primme.numProcs = number_of_nodes;
  primme.procID = this_node;
  primme.globalSumDouble = par_GlobalSumDouble;   // Wrapped function
  primme.matrixMatvec = av_ov;                    // Mat-vec wrapper

  primme_set_method(DEFAULT_MIN_MATVECS, &primme);
  primme.printLevel = 1;
  primme.eps = eig_tol;                   // Maximum residual
  primme.numEvals = Nvec;
  primme.initSize = 0;                    // Number of initial guesses
  if (flag == 1)
    primme.target = primme_smallest;
  else if (flag == -1)
    primme.target = primme_largest;
  else {
    node0_printf("make_evs: Unrecognized flag %d\n", flag);
    terminate(1);
  }
//  primme_display_params(primme);

  // Call the actual EV finder and check return value
  ret = zprimme(eigVal, workVecs, rnorms, &primme);
  while (ret != 0) {
    // Try again with looser residual
    primme.eps *= 10;
    node0_printf("Loosening stopping condition to %.4g\n", primme.eps);
    ret = zprimme(eigVal, workVecs, rnorms, &primme);
  }

  // Copy double-precision temporary fields back into output
  // Each Twist_Fermion has Ndat = 16DIMF non-trivial complex components
  for (ivec = 0; ivec < Nvec; ivec++) {
    iter = Ndat * ivec * sites_on_node;   // This vector in workvecs
    FORALLSITES(i, s) {
      for (j = 0; j < NCOL; j++) {
        for (k = 0; k < NCOL; k++) {
          eigVec[ivec][i].Fsite.e[j][k].real = workVecs[iter].r;
          eigVec[ivec][i].Fsite.e[j][k].imag = workVecs[iter].i;
          iter++;
          FORALLDIR(mu) {
            eigVec[ivec][i].Flink[mu].e[j][k].real = workVecs[iter].r;
            eigVec[ivec][i].Flink[mu].e[j][k].imag = workVecs[iter].i;
            iter++;
          }
          for (mu = 0; mu < NPLAQ; mu++) {
            eigVec[ivec][i].Fplaq[mu].e[j][k].real = workVecs[iter].r;
            eigVec[ivec][i].Fplaq[mu].e[j][k].imag = workVecs[iter].i;
            iter++;
          }
        }
      }
    }
  }
#ifdef DEBUG_CHECK
  if (iter != Nvec * maxn)
    printf("make_evs: iter = %d after output\n", iter);
#endif

  // Print results and check |D^dag D phi - lambda phi|^2
  for (ivec = 0; ivec < Nvec; ivec++) {
    check = 0.0;
    DSq(eigVec[ivec], tmpTF);
    FORALLSITES(i, s) {
      // tTF = tmpTF - eigVal[ivec] * eigVec[ivec]
      scalar_mult_add_TF(&(tmpTF[i]), &(eigVec[ivec][i]),
                                      -1.0 * eigVal[ivec], &tTF);
      check += magsq_TF(&tTF);
    }
    g_doublesum(&check);    // Accumulate across all nodes
    if (flag == 1)  {       // Braces suppress compiler warning
      node0_printf("EIGENVALUE %d %.8g %.8g\n", ivec, eigVal[ivec], check);
    }
    else if (flag == -1)
      node0_printf("BIGEIGVAL  %d %.8g %.8g\n", ivec, eigVal[ivec], check);
  }
  fflush(stdout);

  // Clean up
  free(workVecs);
  free(rnorms);
  primme_Free(&primme);
  return primme.stats.numOuterIterations;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Check matrix elements and eigenvalues of <psi_j | D | psi_i>
// where the psi are eigenvectors of DDag.D
// Have checked that <psi_j | Ddag | psi_i> produces conjugate eigenvalues
void check_Dmat(int Nvec, Twist_Fermion **eigVec) {
  register int i;
  register site *s;
  char N = 'N';
  int ivec, jvec, stat = 0, unit = 1, doub = 2 * Nvec;
  double *store, *work, *eigs, *dum;
  double_complex tc, check;
  Twist_Fermion *tmpTF = malloc(sites_on_node * sizeof(*tmpTF));

  // Allocate double arrays expected by LAPACK
  store = malloc(2 * Nvec * Nvec * sizeof(*store));
  work = malloc(4 * Nvec * sizeof(*work));
  eigs = malloc(2 * Nvec * sizeof(*eigs));
  dum = malloc(2 * sizeof(*dum));

  // Hit each eigVec with D, then contract with every other eigVec
  for (ivec = 0; ivec < Nvec; ivec++) {
    fermion_op(eigVec[ivec], tmpTF, PLUS);

    for (jvec = 0; jvec < Nvec; jvec++) {
      check = cmplx(0.0, 0.0);
      FORALLSITES(i, s) {
        tc = TF_dot(&(eigVec[jvec][i]), &(tmpTF[i]));
        CSUM(check, tc);
      }
      g_dcomplexsum(&check);    // Accumulate across all nodes
//      node0_printf("D[%d, %d] (%.8g, %.4g)\n",
//                   ivec, jvec, check.real, check.imag);

      // Save in column-major double array expected by LAPACK
      store[2 * (jvec + Nvec * ivec)] = check.real;
      store[2 * (jvec + Nvec * ivec) + 1] = check.imag;
    }
  }
  free(tmpTF);

  // Diagonalize <psi|D|psi> using LAPACK
  // Arguments summarized in susy_includes.h
  node0_printf("Using LAPACK to diagonalize <psi_j | D | psi_i>\n");
  zgeev_(&N, &N, &Nvec, store, &Nvec, eigs,
         dum, &unit, dum, &unit, work, &doub, work, &stat);

  // Print resulting eigenvalues
  for (ivec = 0; ivec < Nvec; ivec++)
    node0_printf("D_eig %d (%.6g, %.6g)\n",
                 ivec, eigs[2 * ivec], eigs[2 * ivec + 1]);

  // Free double arrays expected by LAPACK
  free(store);
  free(work);
  free(eigs);
  free(dum);
}
// -----------------------------------------------------------------
