// -----------------------------------------------------------------
// Supersymmetric setup
#include "corr_includes.h"

#define IF_OK if(status==0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and seed, and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
    // stringification kludge from GNU preprocessor manual
    // http://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define XSTR(s) STR(s)
#define STR(s) #s
    // end kludge
    printf("N=4 SYM, Nc = %d, DIMF = %d, fermion rep = adjoint\n",
           NCOL, DIMF);
    printf("Konishi and SUGRA correlator analyses\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");
    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);

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

  this_node = mynode();
  number_of_nodes = numnodes();
  one_ov_N = 1.0 / (Real)NCOL;
  volume = nx * ny * nz * nt;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for fields
void make_fields() {
  int j;
  Real size = 0.0;

  // Temporary matrices and Twist_Fermion
  size += (Real)(3.0 * sizeof(matrix));
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);
  FIELD_ALLOC(staple, matrix);

  // Scalar fields and their bilinear traces
  // Just need one set of N_B / N_K for current measurement
  size += (Real)(N_B * NUMLINK * sizeof(matrix));
  size += (Real)(N_K * NUMLINK * NUMLINK * sizeof(Real));
  for (j = 0; j < N_B; j++)
    FIELD_ALLOC_VEC(Ba[j], matrix, NUMLINK);
  for (j = 0; j < N_K; j++)
    FIELD_ALLOC_MAT(traceBB[j], Real, NUMLINK, NUMLINK);

  // For shifting
  size += (Real)(2.0 * sizeof(Kops));
  FIELD_ALLOC(tempops, Kops);
  FIELD_ALLOC(tempops2, Kops);

#ifdef SMEAR
  // Stout smearing stuff
  size += (Real)(NUMLINK * sizeof(anti_hermitmat));
  FIELD_ALLOC_VEC(Q, anti_hermitmat, NUMLINK);    // To be exponentiated
#endif

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume and seed
  prompt = initial_set();
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
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for Monte Carlo
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status, j;
  Real size;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

#ifdef SMEAR
    // Smearing stuff -- passed to either APE or stout routines by application
    IF_OK status += ask_smear_type(stdin, prompt, &par_buf.smearflag);
    IF_OK status += get_i(stdin, prompt, "Nsmear", &par_buf.Nsmear);
    IF_OK status += get_f(stdin, prompt, "alpha", &par_buf.alpha);
#endif

    // Find out how many lattices we will analyze
    IF_OK status += get_i(stdin, prompt, "Nblock", &par_buf.Nblock);
    IF_OK status += get_i(stdin, prompt, "Nmeas", &par_buf.Nmeas);
    if (par_buf.Nblock * par_buf.Nmeas > MAX_CFG) {
      node0_printf("Error: Can only handle up to %d lattices\n", MAX_CFG);
      node0_printf("       Recompile with different MAX_CFG for more\n");
      status++;
    }

    // Read in lattice names; will reload all in serial
    for (j = 0; j < par_buf.Nblock * par_buf.Nmeas; j++)
      IF_OK status += get_s(stdin, prompt, "reload_serial", par_buf.cfg[j]);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  Nblock = par_buf.Nblock;
  Nmeas = par_buf.Nmeas;
  tot_meas = Nblock * Nmeas;
  for (j = 0; j < tot_meas; j++)
    strcpy(cfg[j], par_buf.cfg[j]);

  // Allocate Konishi and SUGRA operators now that we know Nblock
  ops = malloc(Nblock * sizeof(**ops));
  for (j = 0; j < Nblock; j++)
    FIELD_ALLOC(ops[j], Kops);
  size = (Real)(Nblock * sizeof(Kops)) * sites_on_node;
  node0_printf("\nMallocing %.1f MBytes per core for operators\n", size / 1e6);

#ifdef SMEAR
  smearflag = par_buf.smearflag;
  Nsmear = par_buf.Nsmear;
  alpha = par_buf.alpha;
  if (smearflag == NO_SMEAR) {
    Nsmear = 0;
    alpha = 0.0;
  }
#endif

  // These two arrays only need to be of size total_r <= MAX_pts
  // but allocate MAX_pts memory for simplicity
  // Initialize to nonsense; can be used to check that they have been reset
  lookup = malloc(MAX_pts * sizeof(*lookup));
  norm = malloc(MAX_pts * sizeof(*norm));
  for (j = 0; j < MAX_pts; j++) {
    lookup[j] = -1.0;
    norm[j] = -1.0;
  }

  // Allocate arrays to be used by LAPACK in unit.c
  ipiv = malloc(NCOL * sizeof(*ipiv));
  store = malloc(2 * NCOL * NCOL * sizeof(*store));
  work = malloc(4 * NCOL * sizeof(*work));
  Rwork = malloc((3 * NCOL - 2) * sizeof(*Rwork));
  eigs = malloc(NCOL * sizeof(*eigs));

  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wanted by ../generic/io_lat_utils.c
void write_appl_gauge_info(FILE *fp) {
}
// -----------------------------------------------------------------
