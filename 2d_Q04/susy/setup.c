// -----------------------------------------------------------------
// N=(2,2) SYM setup
#include "susy_includes.h"

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and seed, and send to others
int initial_set() {
  int prompt = 0, status = 0, dir;
  if (mynode() == 0) {
    // Print banner
    printf("N=(2,2) SYM, Nc = %d, DIMF = %d, fermion rep = adjoint\n",
           NCOL, DIMF);
    printf("Microcanonical simulation with refreshing\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else   // Quit!
    printf("Only works for phi algorithm\n");
    exit(1);
#endif
    time_stamp("start");
    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "PBC", &par_buf.PBC);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    // Number of Nth roots to take (in addition to 1/4 power)
    IF_OK status += get_i(stdin, prompt, "Nroot", &par_buf.Nroot);

    // RHMC degree
    IF_OK status += get_i(stdin, prompt, "Norder", &par_buf.Norder);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nx = par_buf.nx;
  nt = par_buf.nt;
  PBC = par_buf.PBC;
  iseed = par_buf.iseed;

  // Lattice volume sanity checks, including dimensional reduction
  length[0] = nx;
  length[1] = nt;
  if (mynode() == 0) {
    FORALLUPDIR(dir) {
      if (length[dir] < 1) {
        printf("nx and nt must both be positive\n");
        exit(1);
      }
    }
  }

  // Set up stuff for RHMC and multi-mass CG
  Nroot = par_buf.Nroot;
  fnorm = malloc(sizeof fnorm * Nroot);
  max_ff = malloc(sizeof max_ff * Nroot);

  Norder = par_buf.Norder;
  amp = malloc(sizeof amp * Norder);
  amp4 = malloc(sizeof amp4 * Norder);
  amp8 = malloc(sizeof amp8 * Norder);
  shift = malloc(sizeof shift * Norder);
  shift4 = malloc(sizeof shift4 * Norder);
  shift8 = malloc(sizeof shift8 * Norder);

  this_node = mynode();
  number_of_nodes = numnodes();
  one_ov_N = 1.0 / (Real)NCOL;
  volume = nx * nt;
  total_iters = 0;
  minus1 = cmplx(-1.0, 0.0);
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for fields
void make_fields() {
#ifdef EIG_POT
  node0_printf("Single-trace scalar potential\n");
#else
  node0_printf("Double-trace scalar potential\n");
#endif
  int Noffdiag = NUMLINK * (NUMLINK - 1);
  Real size = (Real)(2.0 * sizeof(complex));
  FIELD_ALLOC(tr_eta, complex);
  FIELD_ALLOC(tr_dest, complex);

  size += (Real)(2.0 * (2.0 + NUMLINK)) * sizeof(matrix);
  FIELD_ALLOC(site_src, matrix);
  FIELD_ALLOC(site_dest, matrix);
  FIELD_ALLOC_VEC(link_src, matrix, NUMLINK);
  FIELD_ALLOC_VEC(link_dest, matrix, NUMLINK);
  FIELD_ALLOC(plaq_src, matrix);
  FIELD_ALLOC(plaq_dest, matrix);

  // For convenience in calculating action and force
  size += (Real)(2.0 + 3.0 * NUMLINK) * sizeof(matrix);
  size += (Real)(NUMLINK + 3.0 * Noffdiag) * sizeof(complex);
  FIELD_ALLOC(DmuUmu, matrix);
  FIELD_ALLOC(Fmunu, matrix);
  FIELD_ALLOC_VEC(Uinv, matrix, NUMLINK);
  FIELD_ALLOC_VEC(Udag_inv, matrix, NUMLINK);
  FIELD_ALLOC_VEC(UpsiU, matrix, NUMLINK);
  FIELD_ALLOC_VEC(Tr_Uinv, complex, NUMLINK);
  FIELD_ALLOC_MAT_OFFDIAG(plaqdet, complex, NUMLINK);
  FIELD_ALLOC_MAT_OFFDIAG(tempdet, complex, NUMLINK);
  FIELD_ALLOC_MAT_OFFDIAG(ZWstar, complex, NUMLINK);

  // CG Twist_Fermions
  size += (Real)(3.0 * sizeof(Twist_Fermion));
  FIELD_ALLOC(mpm, Twist_Fermion);
  FIELD_ALLOC(pm0, Twist_Fermion);
  FIELD_ALLOC(rm, Twist_Fermion);

  // Temporary matrices and Twist_Fermion
  size += (Real)(3.0 * sizeof(matrix));
  size += (Real)(sizeof(Twist_Fermion));
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);
  FIELD_ALLOC(staple, matrix);
  FIELD_ALLOC(tempTF, Twist_Fermion);

#ifdef CORR
  int j;
  size += (Real)(N_B * NUMLINK * sizeof(matrix));
  size += (Real)(N_K * NUMLINK * NUMLINK * sizeof(Real));
  size += (Real)(2.0 * sizeof(Kops));
  for (j = 0; j < N_B; j++)
    FIELD_ALLOC_VEC(Ba[j], matrix, NUMLINK);
  for (j = 0; j < N_K; j++)
    FIELD_ALLOC_MAT(traceBB[j], Real, NUMLINK, NUMLINK);
  FIELD_ALLOC(tempops, Kops);
  FIELD_ALLOC(tempops2, Kops);
#endif

#ifdef SMEAR
  // Stout smearing stuff
  size += (Real)(NUMLINK * sizeof(anti_hermitmat));
  FIELD_ALLOC_VEC(Q, anti_hermitmat, NUMLINK);    // To be exponentiated
#endif

#if defined(EIG) || defined(PHASE)
  size += (Real)(2.0 * sizeof(Twist_Fermion));
  FIELD_ALLOC(src, Twist_Fermion);
  FIELD_ALLOC(res, Twist_Fermion);
#endif

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
#ifdef PHASE
  // Total number of matvecs is (volume * 4 * DIMF)^2 / 4
  Nmatvecs = volume * 4 * DIMF * volume * DIMF;

  // Total size of matrix is (volume * 4 * DIMF) x (sites_on_node * 4 * DIMF)
  size = (Real)(volume * 4.0 * DIMF * 4.0 * DIMF * sizeof(complex));
  size *= sites_on_node;
  node0_printf("Q has %d columns --> %li matvecs and %.1f MBytes per core...",
               volume * 4 * DIMF, Nmatvecs, size / 1e6);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume and seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Set up offset tables for the five paths
  setup_offset();
  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef SMEAR
// Find out what smearing to use
int ask_smear_type(FILE *fp, int prompt, int *flag) {
  int status = 0;
  char savebuf[256];

  if (prompt != 0)
    printf("enter 'no_smear', 'stout_smear' or 'APE_smear'\n");
  status = fscanf(fp, "%s", savebuf);
  if (status == EOF) {
    printf("ask_smear_type: EOF on STDIN\n");
    return 1;
  }
  if (status != 1) {
    printf("\nask_smear_type: ERROR IN INPUT: ");
    printf("can't read smearing option\n");
    return 1;
  }

  printf("%s\n", savebuf);
  if (strcmp("no_smear", savebuf) == 0)
    *flag = NO_SMEAR;
  else if (strcmp("stout_smear", savebuf) == 0)
    *flag = STOUT_SMEAR;
  else if (strcmp("APE_smear", savebuf) == 0)
    *flag = APE_SMEAR;
  else {
    printf("Error in input: invalid smear type\n");
    printf("Only no_smear, stout_smear and APE_smear supported\n");
    return 1;
  }
  return 0;
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for Monte Carlo
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;
#ifdef EIG
  int i;
#endif
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);

    // Number of fermion and gauge steps
    IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps[0]);
    IF_OK status += get_i(stdin, prompt, "nstep_gauge", &par_buf.nsteps[1]);

    // Trajectories between propagator measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.propinterval);

    // lambda, kappa_u1, bmass, fmass, G
    IF_OK status += get_f(stdin, prompt, "lambda", &par_buf.lambda);
    IF_OK status += get_f(stdin, prompt, "kappa_u1", &par_buf.kappa_u1);
    IF_OK status += get_f(stdin, prompt, "bmass", &par_buf.bmass);
    IF_OK status += get_f(stdin, prompt, "fmass", &par_buf.fmass);
    IF_OK status += get_f(stdin, prompt, "G", &par_buf.G);

#ifdef SMEAR
    // Smearing stuff -- passed to either APE or stout routines by application
    IF_OK status += ask_smear_type(stdin, prompt, &par_buf.smearflag);
    IF_OK status += get_i(stdin, prompt, "Nsmear", &par_buf.Nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha", &par_buf.alpha);
#endif

    // Maximum conjugate gradient iterations
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);

    // Error per site for conjugate gradient
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x * x;
    }

#ifdef BILIN
    // Number of stochastic sources for fermion bilinear and susy trans
    // Also used for stochastic mode number computation
    IF_OK status += get_i(stdin, prompt, "nsrc", &par_buf.nsrc);
#endif

#ifdef EIG
    // Number of eigenvalues to calculate
    IF_OK status += get_i(stdin, prompt, "Nvec", &par_buf.Nvec);
    IF_OK status += get_f(stdin, prompt, "eig_tol", &par_buf.eig_tol);
    IF_OK status += get_i(stdin, prompt, "maxIter", &par_buf.maxIter);
#endif

#ifdef PHASE
    // Optional checkpointing for pfaffian computation
    IF_OK status += get_i(stdin, prompt, "ckpt_load", &par_buf.ckpt_load);
    IF_OK status += get_i(stdin, prompt, "ckpt_save", &par_buf.ckpt_save);
#endif

#ifdef WLOOP
    // Find out whether or not to gauge fix to Coulomb gauge
    IF_OK status += ask_gauge_fix(stdin, prompt, &par_buf.fixflag);
#endif

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
                                       par_buf.savefile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length = par_buf.traj_length;
  nsteps[0] = par_buf.nsteps[0];
  nsteps[1] = par_buf.nsteps[1];

  propinterval = par_buf.propinterval;
  fixflag = par_buf.fixflag;
  niter = par_buf.niter;
  rsqmin = par_buf.rsqmin;

  lambda = par_buf.lambda;
  kappa_u1 = par_buf.kappa_u1;
  bmass = par_buf.bmass;
  fmass = par_buf.fmass;
  G = par_buf.G;
  if (G > IMAG_TOL)
    doG = 1;
  else
    doG = 0;

  kappa = (Real)NCOL * 0.25 / lambda;
  node0_printf("lambda=%.4g --> kappa=Nc/(4lambda)=%.4g\n",
               lambda, kappa);
  node0_printf("C2=%.4g\n", C2);    // Currently hardwired in defines.h

#ifdef SMEAR
  smearflag = par_buf.smearflag;
  Nsmear = par_buf.Nsmear;
  alpha = par_buf.alpha;
  if (smearflag == NO_SMEAR) {
    Nsmear = 0;
    alpha = 0.0;
  }
#endif

#ifdef BILIN
  nsrc = par_buf.nsrc;
#endif

#ifdef EIG
  // Include some mallocs here (make_fields has already been called)
  // (This memory usage will not be reported to the user)
  Nvec = par_buf.Nvec;
  eigVal = malloc(sizeof *eigVal * Nvec);
  eigVec = malloc(sizeof *eigVec * Nvec);
  for (i = 0; i < Nvec; i++)
    FIELD_ALLOC(eigVec[i], Twist_Fermion);

  eig_tol = par_buf.eig_tol;
  maxIter = par_buf.maxIter;
#endif

#ifdef PHASE
  ckpt_load = par_buf.ckpt_load;
  ckpt_save = par_buf.ckpt_save;
#endif

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);

  // Allocate arrays to be used by LAPACK in determinant.c
  // Needs to be above compute_Uinv()
  ipiv = malloc(sizeof *ipiv * NCOL);
  store = malloc(sizeof *store * 2 * NCOL * NCOL);
  work = malloc(sizeof *work * 4 * NCOL);

  // Allocate some more arrays to be used by LAPACK in unit.c
  Rwork = malloc(sizeof *Rwork * (3 * NCOL - 2));
  eigs = malloc(sizeof *eigs * NCOL);

  // Compute initial plaqdet, DmuUmu and Fmunu
  compute_plaqdet();
  compute_Uinv();
  compute_DmuUmu();
  compute_Fmunu();
  return 0;
}
// -----------------------------------------------------------------
