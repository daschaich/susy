// -----------------------------------------------------------------
// N=4 SYM setup
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
    printf("N=4 SYM, Nc = %d, DIMF = %d, fermion rep = adjoint\n",
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
#ifdef DIMREDUCE
    printf("Dimensionally reduced calculation\n");
#endif
    time_stamp("start");
    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
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
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  PBC = par_buf.PBC;
  iseed = par_buf.iseed;

  // Lattice volume sanity checks, including dimensional reduction
  length[0] = nx;
  length[1] = ny;
  length[2] = nz;
  length[3] = nt;
  if (mynode() == 0) {
    FORALLUPDIR(dir) {
      if (length[dir] < 1) {
        printf("{nx, ny, nz, nt} must all be positive\n");
        exit(1);
      }
    }
#ifdef DIMREDUCE
    if (nx > 1 && ny > 1 && nz > 1 && nt > 1) {
      printf("WARNING: Compiled with dimensional reduction ");
      printf("but running without any reduced dims\n");
    }
#else
    if (nx == 1 || ny == 1 || nz == 1 || nt == 1) {
      printf("WARNING: Running with reduced dim(s) ");
      printf("but didn't compile with -DDIMREDUCE\n");
    }
#endif
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
  volume = nx * ny * nz * nt;
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
#ifdef RESCALE
  node0_printf("Rescaled fermion operator\n");
#endif
  Real size = (Real)(2.0 * sizeof(complex));
  FIELD_ALLOC(tr_eta, complex);
  FIELD_ALLOC(tr_dest, complex);

  size += (Real)(2.0 * (1.0 + NUMLINK + NPLAQ) * sizeof(matrix));
  FIELD_ALLOC(site_src, matrix);
  FIELD_ALLOC(site_dest, matrix);
  FIELD_ALLOC_VEC(link_src, matrix, NUMLINK);
  FIELD_ALLOC_VEC(link_dest, matrix, NUMLINK);
  FIELD_ALLOC_VEC(plaq_src, matrix, NPLAQ);
  FIELD_ALLOC_VEC(plaq_dest, matrix, NPLAQ);

  // For convenience in calculating action and force
  size += (Real)(1.0 + NPLAQ + 3.0 * NUMLINK) * sizeof(matrix);
  size += (Real)(NUMLINK + 6.0 * NPLAQ) * sizeof(complex);
  FIELD_ALLOC(DmuUmu, matrix);
  FIELD_ALLOC_VEC(Fmunu, matrix, NPLAQ);
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
  size += (Real)(sizeof(gather_mat));
  size += (Real)((NUMLINK + 2.0) * sizeof(gather_vec));
  size += (Real)(sizeof(Twist_Fermion));
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);
  FIELD_ALLOC(staple, matrix);
  FIELD_ALLOC(tempgathmat, gather_mat);
  FIELD_ALLOC(tempgathvec, gather_vec);
  FIELD_ALLOC(tempgathvec2, gather_vec);
  FIELD_ALLOC_VEC(tempgathvec3, gather_vec, NUMLINK);
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

#if defined(CHEB) || defined(MODE)
  // For Z2 random source
  size += (Real)(sizeof(Twist_Fermion));
  FIELD_ALLOC(z_rand, Twist_Fermion);
#endif

#ifdef MODE
  // Temporary TF for stochastic mode number
  size += (Real)(5.0 * sizeof(Twist_Fermion));
  FIELD_ALLOC(XPXSq, Twist_Fermion);
  FIELD_ALLOC(hX, Twist_Fermion);
  FIELD_ALLOC(dest, Twist_Fermion);
  FIELD_ALLOC(bj, Twist_Fermion);
  FIELD_ALLOC(bjp1, Twist_Fermion);
#endif

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
#ifdef PHASE
  // Total number of matvecs is (volume * 16 * DIMF)^2 / 4
  Nmatvecs = volume * 16 * DIMF * volume * 4 * DIMF;

  // Total size of matrix is (volume * 16 * DIMF) x (sites_on_node * 16 * DIMF)
  size = (Real)(volume * 16.0 * DIMF * 16.0 * DIMF * sizeof(complex));
  size *= sites_on_node;
  node0_printf("Q has %d columns --> %li matvecs and %.1f MBytes per core...",
               volume * 16 * DIMF, Nmatvecs, size / 1e6);
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
  // Set up offset tables for gathers
  setup_offset();
  setup_qclosed_offset();
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
#if defined(EIG) || defined(MODE)
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
#ifdef DIMREDUCE
    IF_OK status += get_f(stdin, prompt, "cWline", &par_buf.cWline);
#endif

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

#if defined(CHEB) || defined(MODE)
    // Number of stochastic sources
    IF_OK status += get_i(stdin, prompt, "Nstoch", &par_buf.Nstoch);
#endif

#ifdef CHEB
    // How many Chebyshev coefficients to compute
    IF_OK status += get_i(stdin, prompt, "cheb_order", &par_buf.cheb_order);

    // Bounds on spectral range
    IF_OK status += get_f(stdin, prompt, "lambda_min", &par_buf.lambda_min);
    IF_OK status += get_f(stdin, prompt, "lambda_max", &par_buf.lambda_max);
#endif

#ifdef MODE
    // Which order polynomial to use in step function
    IF_OK status += get_i(stdin, prompt, "step_order", &par_buf.step_order);

    // A maximum of MAX_OMEGA points at which to evaluate the mode number
    IF_OK status += get_i(stdin, prompt, "numOmega", &par_buf.numOmega);
    if (par_buf.numOmega > MAX_OMEGA) {
      node0_printf("ERROR: Need to recompile for numOmega > %d\n",
                   MAX_OMEGA);
      status++;
    }
    for (i = 0; i < par_buf.numOmega; i++)
      IF_OK status += get_f(stdin, prompt, "Omega", &par_buf.Omega[i]);
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

#ifdef DIMREDUCE
  cWline = par_buf.cWline;
#endif

  kappa = (Real)NCOL * 0.25 / lambda;
  node0_printf("lambda=%.4g --> kappa=Nc/(4lambda)=%.4g\n",
               lambda, kappa);
  node0_printf("C2=%.4g\n", C2);    // Currently hardwired in defines.h

#ifdef WLOOP
  if (fixflag != NO_GAUGE_FIX)
    FIELD_ALLOC(gfix_u, su2_matrix);
#endif

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

#if defined(CHEB) || defined(MODE)
  Nstoch = par_buf.Nstoch;

  // Normalization factor for errors from averaging over stochastic sources
  if (Nstoch > 1)
    sqrt1_ov_Nm1 = 1.0 / sqrt((Real)Nstoch - 1.0);
  else
    sqrt1_ov_Nm1 = 0.0;
#endif

#ifdef CHEB
  cheb_order = par_buf.cheb_order;
  cheb_coeff = malloc(sizeof(Real) * cheb_order);
  cheb_err = malloc(sizeof(Real) * cheb_order);

  lambda_min = par_buf.lambda_min;
  lambda_max = par_buf.lambda_max;
#endif

#ifdef MODE
  // Save sources to reuse for each Omega
  source = malloc(sizeof *source * Nstoch);
  for (i = 0; i < Nstoch; i++)
    FIELD_ALLOC(source[i], Twist_Fermion);

  step_order = par_buf.step_order;
  step_coeff = malloc(sizeof(Real) * (step_order + 1));

  numOmega = par_buf.numOmega;
  mode = malloc(sizeof(Real) * numOmega);
  err = malloc(sizeof(Real) * numOmega);
  Omega = malloc(sizeof(Real) * numOmega);
  for (i = 0; i < numOmega; i++)
    Omega[i] = par_buf.Omega[i];
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
