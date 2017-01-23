// -----------------------------------------------------------------
// Main procedure for N=4 SYM modenumber calculation
// It includes two different methods for the determination
// one is the standard Giusti Luescher method, the other one is
// based on an expansion in Chebyshev polynomials.
// For the Giusti Luescher method an optimized polynomial for the
// step function has to be provided in coeff.txt.
// The modenumber is in that case calculated for a set of shifts
// specified in shifts.txt.
// The Chebyshev approximation method requires some bound on the
// minimal and maximal eigenvalues. It provides a set of coefficients
// for the Chenyshev expansion of the eigenvalue density. These need to be
// post-processed to calculate the modenumber.
#define CONTROL
#include "susy_includes.h"
// -----------------------------------------------------------------


// Reading of the shift and coefficient parameter files. This has to be done
// on for each process since a broadcast is so far not included.
void readNumberFile(unsigned int* order,Real** coeff,const char* filename){
	unsigned int count=0;
	double tmp;
	FILE *file;
    *order=0;
    (*coeff)=NULL;
	file = fopen(filename, "r");
	if (file == NULL){
		printf("ERROR: File %s not readable\n",filename);
		return;
	}
	while(fscanf(file, "%lg", &tmp)==1){
		count++;
	}
	if(count==0){
		printf("ERROR: File %s is empty\n",filename);
		return;
	}
	file = freopen(filename, "r",file);
	*coeff=malloc(count*sizeof(Real));
	*order=count;
	count=0;
	while(fscanf(file, "%lg", &tmp)==1){

		(*coeff)[count]=tmp;
		count++;
	}
}


// -----------------------------------------------------------------
// Main measurement: as for the eigenvalue, it includes some additional cheap measurements for cross checks.
int main(int argc, char *argv[]) {
  int prompt, dir;
  double ss_plaq, st_plaq, dtime, plpMod = 0.0;
  double linktr[NUMLINK], linktr_ave, linktr_width;
  double link_det[NUMLINK], det_ave, det_width;
  complex plp = cmplx(99.0, 99.0);
	int k;
#ifdef SEPARATE_FILE
	FILE * outfile;
#endif
	Real* ckcoeff = NULL;

	Real* shiftparam=NULL;
	unsigned int noshifts=0;
    Real* coeff=NULL;
    unsigned int coefforder;
    Real* result=NULL;
    Real* resultstaterr=NULL;

#ifndef MODE
#error "Don't use control_mode unless compiling with -DMODE!"
#endif

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  epsilon();
  setup_PtoP();
  setup_FQ();

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial plaquette and bosonic action
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g ", ss_plaq, st_plaq, ss_plaq + st_plaq);
  ss_plaq = gauge_action(NODET);
  node0_printf("%.8g\n", ss_plaq / (double)volume);

  // Do "local" measurements to check configuration
  // Tr[Udag.U] / N
  linktr_ave = link_trace(linktr, &linktr_width,
                          link_det, &det_ave, &det_width);
  node0_printf("FLINK");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", linktr[dir]);
  node0_printf(" %.6g %.6g\n", linktr_ave, linktr_width);
  node0_printf("FLINK_DET");
  for (dir = XUP; dir < NUMLINK; dir++)
    node0_printf(" %.6g", link_det[dir]);
  node0_printf(" %.6g %.6g\n", det_ave, det_width);

  // Polyakov loop and plaquette measurements
  // Format: GMES Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq
  plp = ploop(TUP, NODET, &plpMod);
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("GMES %.8g %.8g 0 %.8g %.8g ",
               plp.real, plp.imag, ss_plaq, st_plaq);

  // Bosonic action (printed twice by request)
  // Might as well spit out volume average of Polyakov loop modulus
  ss_plaq = gauge_action(NODET) / (double)volume;
  node0_printf("%.8g ", ss_plaq);
  node0_printf("%.8g\n", plpMod);
  node0_printf("BACTION %.8g\n", ss_plaq);

  // deacitivate measurement by setting Nest or order of the polynomial to zero.
  if(step_order>0 && Nest>0){
	  // Main calculation in Chebyshev method
     calculate_coeff(Nest, step_order, lambda_min, lambda_max, &ckcoeff);

#ifdef SEPARATE_FILE
	if (this_node == 0) {
		outfile = fopen("chebcoeff.txt", "a");
		for (k = 0; k < step_order; k++) {
			fprintf(outfile, "%d %d %.16g\n", 0, k, ckcoeff[k]);
		}
		fclose(outfile);
	}
#else
	// printout
	for (k = 0; k < step_order; k++) {
		node0_printf("MODE_CHEBYSHEV_CKCOEFF: %d %d %.16lg\n", 0, k, ckcoeff[k]);
	}
#endif
  }else{
	  node0_printf("Chebyshev coefficient calculation deactivated.\n");
  }
    node0_printf("Reading files \n");
	readNumberFile(&noshifts,&shiftparam,"shifts.txt");
	readNumberFile(&coefforder,&coeff,"coeff.txt");

  // Measurement of Giusti Luescher method, deactivate if coefficient or shift file not present.
  if(noshifts==0 || coefforder==0){
	printf("ERROR: no coefficients or shifts avaiable.\n");
  }else{
	  node0_printf("Shifts: ");
	for (k = 0; k < noshifts; k++) {
		node0_printf("%.8lg ",shiftparam[k]);
	}
	  node0_printf("\n");
	  node0_printf("Coeff: ");
	for (k = 0; k < coefforder; k++) {
		node0_printf("%.8lg ",coeff[k]);
	}
	  node0_printf("\n");
	result=malloc(noshifts*sizeof(Real));
	resultstaterr=malloc(noshifts*sizeof(Real));
	// Main calculation GL method.
	calculateGLMethod(noshifts,shiftparam,rescalefact_gl,result,resultstaterr,nest_gl,coefforder,coeff,epsilon_gl,residgoal_gl,maxit_gl);
#ifdef SEPARATE_FILE
	if (this_node == 0) {
		outfile = fopen("chebmodegl.txt", "a");
		for (k = 0; k < noshifts; k++) {
			fprintf(outfile, "%d %d %.16g %.16g %.16g\n", 0, k,shiftparam[k], result[k],resultstaterr[k]);
		}
		fclose(outfile);
	}
#else
	// printout
	for (k = 0; k < noshifts; k++) {
		node0_printf("MODE_GL_MODENUM: %d %d %.16lg %.16lg %.16lg\n", 0, k,shiftparam[k], result[k],resultstaterr[k]);
	}
#endif
  }

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n", total_iters);
  fflush(stdout);
	free(ckcoeff);
  if(shiftparam){
   free(shiftparam);
  }
  if(coeff){
    free(coeff);
  }
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
