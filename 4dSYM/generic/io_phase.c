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
void loadQ(complex **Q, complex *diag, int ckpt_load) {
  int i, j, Ndat = 16 * DIMF;
  char infile[80];
  FILE *fp = NULL, *fd = NULL;

  node0_printf("Reloading columns %d--%d\n", ckpt_load + 1, volume * Ndat);

  // First load remaining columns of Q
  sprintf(infile, "%s.Q%d-node%d", startfile, ckpt_load, this_node);
  fp = fopen(infile, "r");    // Open to read
  if (fp == NULL) {
    printf("loadQ: node%d can't open file %s\n", this_node, infile);
    fflush(stdout);
    terminate(1);
  }

  for (i = ckpt_load; i < volume * Ndat; i++) {
    for (j = 0; j < sites_on_node * Ndat; j++)
      fscanf(fp, "%lg\n%lg\n", &(Q[i][j].real), &(Q[i][j].imag));
  }
  fclose(fp);

  // Now load completed diagonal elements only on node0
  if (this_node == 0) {
    sprintf(infile, "%s.diag%d", startfile, ckpt_load);
    fd = fopen(infile, "r");    // Open to read
    if (fd == NULL) {
      printf("saveQ: node0 can't open file %s\n", infile);
      fflush(stdout);
      terminate(1);
    }

    for (i = 1; i < ckpt_load; i += 2)  // Every other is trivial
      fscanf(fd, "%lg\n%lg\n", &(diag[i].real), &(diag[i].imag));
    fclose(fd);
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
void saveQ(complex **Q, complex *diag, int ckpt_save) {
  int i, j, Ndat = 16 * DIMF;
  char outfile[80];
  FILE *fp = NULL, *fd = NULL;

  node0_printf("Dumping columns %d--%d\n", ckpt_save + 1, volume * Ndat);

  // First save completed diagonal elements only on node0
  if (this_node == 0) {
    sprintf(outfile, "%s.diag%d", startfile, ckpt_save);
    fd = fopen(outfile, "w");    // Open to write
    if (fd == NULL) {
      printf("saveQ: node%d can't open file %s\n", this_node, outfile);
      fflush(stdout);
      terminate(1);
    }

    for (i = 1; i < ckpt_save; i += 2)  // Every other is trivial
      fprintf(fd, "%.16g\n%.16g\n", diag[i].real, diag[i].imag);
    fclose(fd);
  }

  // Then save remaining columns of Q
  sprintf(outfile, "%s.Q%d-node%d", startfile, ckpt_save, this_node);
  fp = fopen(outfile, "w");    // Open to write
  if (fp == NULL) {
    printf("saveQ: node%d can't open file %s\n", this_node, outfile);
    fflush(stdout);
    terminate(1);
  }

  // This will save lots of zeros on some nodes
  // but it's the easy way to keep from missing data on other nodes
  for (i = ckpt_save; i < volume * Ndat; i++) {
    for (j = 0; j < sites_on_node * Ndat; j++)
      fprintf(fp, "%.16g\n%.16g\n", Q[i][j].real, Q[i][j].imag);
    free(Q[i]);
  }
  fclose(fp);
}
#endif
// -----------------------------------------------------------------
