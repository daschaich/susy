// -----------------------------------------------------------------
// Main procedure for BFSS/BMN evolution and measurements
#define CONTROL
#include "susy_includes.h"

int main(int argc, char *argv[]) {
  int prompt, j, h, l, k, i;
  int traj_done, s_iters, avs_iters = 0, avm_iters = 0, Nmeas = 0;
  Real f_eps, b_eps;
  double b_act, dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
  double ave_eigs[NCOL], eig_widths[NCOL], min_eigs[NCOL], max_eigs[NCOL];
  complex plp = cmplx(99.0, 99.0);
  double re, im;
  register site *s;

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  setup_gamma();
  setup_rhmc();
    
  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();
    
  // The part below is added for testing with the config generated from serial C++ code
  // We will use the same "fresh" option and then overwrite the link, scalars at this point
  // Now only for N = 2 . Should suffice for testing with serial C++ code for now.
    
  #ifndef PUREGAUGE
  node0_printf("ERROR in reading config because NO FERMIONS, aborting\n");
  terminate(1);
  #endif
    
  double** data=malloc(1000*sizeof(double*));
  for(i=0;i<1000;i++)
  data[i]=malloc(1000*sizeof(double));
    
  FILE *input;
  input = fopen("config.txt","r");
    
  if (input == NULL) {
      fprintf(stderr, "Can't open file %s!\n","CONFIG");
      exit(1);
  }
    
  for(i = 0; i < NCOL*NCOL*2; i++) {
     for(j = 0; j < 10*nt; j++){
            if (!fscanf(input, "%lf", &data[i][j]))
            break;
      }
  }
  fclose(input);
    
  FORALLSITES(i, s){
        for (j = 0; j < NCOL; j++) {
        for (k = 0; k < NCOL; k++) {
            
        re = data[0][4*j + 2*k + i*NCOL*NCOL*2];
        im = data[0][4*j + 2*k + 1 + i*NCOL*NCOL*2];
        //printf(" Site is [%i] and Data [%i][%i] = %lf\n", i, 0, 4*j + 2*k + i*8, re);
        s->link.e[j][k] = cmplx(re, im);
        }
        }
  }
    
  FORALLSITES(i, s) {
    for (l = 1; l < NSCALAR + 1 ; l++) {
        clear_mat(&(s->X[l]));
        for (j = 0; j < NCOL; j++) {
            for (k = 0; k < NCOL; k++) {
            re = data[l][4*j + 2*k + i*NCOL*NCOL*2];
            im = data[l][4*j + 2*k + 1 + i*NCOL*NCOL*2];
            s->X[l].e[j][k] = cmplx(re, im);
            }
        }
    }
}
    
node0_printf("CONFIG. FROM SERIAL CODE LOADED OVERWRITING ALL BEFORE \n");

  // Check: compute initial bosonic action
  b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]));
  node0_printf("START %.8g\n", b_act / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  // Perform warmup trajectories
  f_eps = traj_length / (Real)nsteps[0];
  b_eps = f_eps / (Real)(2 * nsteps[1]);
  node0_printf("f_eps %.4g b_eps %.4g\n", f_eps, b_eps);
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories, reunitarizations and measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    s_iters = update();
    avs_iters += s_iters;

    // Do "local" measurements every trajectory!
    // Tr[Udag.U] / N
    Xtr_ave = scalar_trace(Xtr, &Xtr_width);
    node0_printf("SCALAR SQUARES");
    for (j = 0; j < NSCALAR; j++)
      node0_printf(" %.6g", Xtr[j]);

    node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);
    // Polyakov loop measurement
    // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters
    plp = ploop();
    node0_printf("GMES %.8g %.8g %d ", plp.real, plp.imag, s_iters);

    // Bosonic action
    b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]));
    node0_printf("%.8g %.8g %.8g %.8g\n",
                 b_act / (double)nt, Xtr[0] / (double)nt,
                 Xtr[1] / (double)nt, Xtr[2] / (double)nt);

    // Monitor scalar eigenvalues
    // Format: SCALAR_EIG # ave width min max
    scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
    for (j = 0; j < NCOL; j++) {
      node0_printf("SCALAR_EIG %d %.6g %.6g %.6g %.6g\n",
                   j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
    }

    // Less frequent measurements every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
#ifdef CORR
      // Konishi and SUGRA
      konishi();
#endif

#ifdef BILIN
      // Ward identity involving eta.psi_a fermion bilinear
      Nmeas++;
      avm_iters += bilinearWard();
#endif
    }
    fflush(stdout);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final bosonic action
  b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]));
  node0_printf("STOP %.8g\n", b_act / (double)nt);

  node0_printf("Average CG iters for steps: %.4g\n",
               (double)avs_iters / trajecs);
  if (Nmeas > 0) {
    node0_printf("Average CG iters for measurements: %.4g\n",
                 (double)avm_iters / Nmeas);
  }
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  // total_iters is accumulated in the multiCG itself
  // Should equal total for steps plus measurements
  node0_printf("total_iters = %d\n\n", total_iters);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
