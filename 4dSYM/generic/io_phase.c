// -----------------------------------------------------------------
// Crude IO for diagonal elements and columns of Q in pfaffian computation
// For now keep things simple by not bothering with a header or checksums
// In fact, let's not even bother with binary reading/writing for now
// Each node reads and writes its own ascii file
#include "generic_includes.h"
#include "../include/io_lat.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void load_diag(complex *diag, int ckpt_load) {
  int i;
  char infile[80];
  FILE *fp = NULL;

  // Only load completed diagonal elements on node0
  if (this_node == 0) {
    sprintf(infile, "%s.diag%d", startfile, ckpt_load);
    fp = fopen(infile, "r");    // Open to read
    if (fp == NULL) {
      printf("saveQ: node0 can't open file %s\n", infile);
      fflush(stdout);
      terminate(1);
    }

    for (i = 1; i < ckpt_load; i += 2)  // Every other is trivial
      fscanf(fp, "%lg\n%lg\n", &(diag[i].real), &(diag[i].imag));
    fclose(fp);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void loadQ(complex **Q, int ckpt_load) {
  int i, j, ret = 1, Ndat = 16 * DIMF;
  char infile[80];
  Real re, im;
  FILE *fp = NULL;

  node0_printf("Reloading columns %d--%d\n", ckpt_load + 1, volume * Ndat);
  sprintf(infile, "%s.Q%d-node%d", startfile, ckpt_load, this_node);
  fp = fopen(infile, "r");    // Open to read
  if (fp == NULL) {
    printf("loadQ: node%d can't open file %s\n", this_node, infile);
    fflush(stdout);
    terminate(1);
  }

  i = 1;
  while (ret != EOF && i > 0) {
    ret = fscanf(fp, "%d\t%d\t%lg\t%lg\n", &i, &j, &re, &im);
    if (ret == EOF || i == -1)    // End of file key
      break;
    Q[i][j].real = re;
    Q[i][j].imag = im;
  }
  // Check that we reached the end
  // If not, that's bad since the previous round will have to be re-run
  // TODO: Replace this with an actual checksum
  if (i != -1) {
    printf("loadQ: node%d didn't reach the end of Q\n", this_node);
    fflush(stdout);
    terminate(1);
  }
  fclose(fp);
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void save_diag(complex *diag, int ckpt_save) {
  int i;
  char outfile[80];
  FILE *fp = NULL;

  // Only save completed diagonal elements on node0
  if (this_node == 0) {
    sprintf(outfile, "%s.diag%d", startfile, ckpt_save);
    fp = fopen(outfile, "w");    // Open to write
    if (fp == NULL) {
      printf("saveQ: node%d can't open file %s\n", this_node, outfile);
      fflush(stdout);
      terminate(1);
    }

    for (i = 1; i < ckpt_save; i += 2)  // Every other is trivial
      fprintf(fp, "%.16g\n%.16g\n", diag[i].real, diag[i].imag);
    fclose(fp);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void saveQ(complex **Q, int ckpt_save) {
  int i, j, Ndat = 16 * DIMF;
  char outfile[80];
  FILE *fp = NULL;

  node0_printf("Dumping columns %d--%d\n", ckpt_save + 1, volume * Ndat);
  sprintf(outfile, "%s.Q%d-node%d", startfile, ckpt_save, this_node);
  fp = fopen(outfile, "w");    // Open to write
  if (fp == NULL) {
    printf("saveQ: node%d can't open file %s\n", this_node, outfile);
    fflush(stdout);
    terminate(1);
  }

  // Try to avoid reading and writing lots of zeros
  for (i = ckpt_save; i < volume * Ndat; i++) {
    for (j = 0; j < sites_on_node * Ndat; j++)
      if (Q[i][j].real != 0.0 || Q[i][j].real != 0.0) { // Loose condition
        fprintf(fp, "%d\t%d\t%.16g\t%.16g\n",
                    i, j, Q[i][j].real, Q[i][j].imag);
    }
    free(Q[i]);
  }
  fprintf(fp, "-1\t-1\t-1\t-1\n");  // To check that we reached the end
  fclose(fp);
}
#endif
// -----------------------------------------------------------------

