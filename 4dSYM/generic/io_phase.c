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
// Only load completed diagonal elements on node0
void load_diag(complex *diag, int ckpt_load) {
  int i;
  char infile[80];
  FILE *fp = NULL;

  if (this_node == 0) {
    sprintf(infile, "%s.diag%d", startfile, ckpt_load);
    fp = fopen(infile, "r");    // Open to read
    if (fp == NULL) {
      printf("load_diag: node0 can't open file %s\n", infile);
      fflush(stdout);
      terminate(1);
    }

    for (i = 1; i < ckpt_load; i += 2)  // Every other is trivial
      fscanf(fp, "%lg\n%lg\n", &(diag[i].real), &(diag[i].imag));
    fclose(fp);
  }
  g_sync();   // Don't let other nodes race ahead
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
// Only load non-zero elements of remaining columns of Q
// Read from node0 only -- transfer data to all other nodes one by one
void loadQ(complex **Q, int ckpt_load) {
  int i, j, savi = 0, savj = 0, ret = 1, Ndat = 16 * DIMF, Qlen;
  int mynode = 0, destnode = 0, Nlines = 0;
  char infile[80];
  Real re, im;
  complex **tbuf = NULL;   // Only malloc on node0 if numnodes() > 1
  FILE *fp = NULL;

  // Total length of each column of Q
  Qlen = sites_on_node * Ndat * sizeof(complex);

  if (this_node == 0) {
    node0_printf("Reloading columns %d--%d\n", ckpt_load + 1, volume * Ndat);
    sprintf(infile, "%s.Q%d", startfile, ckpt_load);
    fp = fopen(infile, "r");    // Open to read
    if (fp == NULL) {
      printf("loadQ: node%d can't open file %s\n", this_node, infile);
      fflush(stdout);
      terminate(1);
    }

    // Allocate temporary storage to be passed to other nodes
    // Only need if numnodes() > 1, but hard-code it in any case
    // to suppress compiler warnings
    tbuf = malloc(volume * Ndat * sizeof(complex*));
    for (i = ckpt_load; i < volume * Ndat; i++) {
      tbuf[i] = malloc(sites_on_node * Ndat * sizeof(complex));
      for (j = 0; j < sites_on_node * Ndat; j++)
        tbuf[i][j] = cmplx(0.0, 0.0);
    }
    if (tbuf[volume * Ndat - 1] == NULL) {
      printf("loadQ: can't malloc tbuf\n");
      fflush(stdout);
      terminate(1);
    }

    // Fill Q[i][j] on node0
    while (ret != EOF && i > 0) {
      ret = fscanf(fp, "%d\t%d\t%d\t%lg\t%lg\n", &mynode, &i, &j, &re, &im);
      Nlines++;
      if (mynode == 1) {
        tbuf[i][j] = cmplx(re, im);   // So we don't lose this
        break;
      }
      else if (mynode == -1) {
        savi = i;
        break;
      }
      Q[i][j] = cmplx(re, im);
    }
  }
  g_sync();   // Don't let other nodes race ahead

  // Now fill in each other node in turn
  for (destnode = 1; destnode < numnodes(); destnode++) {
    // Fill rest of buffer on node0
    if (this_node == 0) {
      while (ret != EOF && i > 0) {
        ret = fscanf(fp, "%d\t%d\t%d\t%lg\t%lg\n", &mynode, &i, &j, &re, &im);
        Nlines++;
        if (mynode != destnode) {
          savi = i;   // May still need to write below
          savj = j;
          break;
        }
        tbuf[i][j] = cmplx(re, im);
      }
    }
    g_sync();   // Don't let other nodes race ahead

    // Now send buffer to destnode
    for (i = ckpt_load; i < volume * Ndat; i++) {
      if (this_node == 0)
        send_field((char *)tbuf[i], Qlen, destnode);
      else if (this_node == destnode)
        get_field((char *)Q[i], Qlen, 0);

      g_sync();   // To be cautious
    }

    // Reset tbuf if we have another node to do
    if (this_node == 0 && mynode != -1) {
      for (i = ckpt_load; i < volume * Ndat; i++) {
        for (j = 0; j < sites_on_node * Ndat; j++)
          tbuf[i][j] = cmplx(0.0, 0.0);
      }

      tbuf[savi][savj] = cmplx(re, im);
    }
    g_sync();   // Don't let other nodes race ahead
  }

  // Now we should be done
  if (this_node == 0) {
    if (mynode != -1) {
      printf("loadQ: should be done but tag = %d\n", mynode);
      fflush(stdout);
      terminate(1);
    }

      // Check that the correct number of elements were read
    if (Nlines - 1 != savi) {   // End of file key
      printf("loadQ: read %d elements but should have had %d\n",
             Nlines - 1, savi);
      fflush(stdout);
      terminate(1);
    }

    fclose(fp);
    for (i = ckpt_load; i < volume * Ndat; i++)
      free(tbuf[i]);
    free(tbuf);
  }
  g_sync();   // Don't let other nodes race ahead
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
// Only save completed diagonal elements on node0
void save_diag(complex *diag, int ckpt_save) {
  int i;
  char outfile[80];
  FILE *fp = NULL;

  if (this_node == 0) {
    sprintf(outfile, "%s.diag%d", startfile, ckpt_save);
    fp = fopen(outfile, "w");    // Open to write
    if (fp == NULL) {
      printf("save_diag: node%d can't open file %s\n", this_node, outfile);
      fflush(stdout);
      terminate(1);
    }

    for (i = 1; i < ckpt_save; i += 2)  // Every other is trivial
      fprintf(fp, "%.16g\n%.16g\n", diag[i].real, diag[i].imag);
    fflush(fp);
    fclose(fp);
  }
  g_sync();   // Don't let other nodes race ahead
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef PHASE
// Only save non-zero elements of remaining columns of Q
// Write from node0 only -- transfer data from all other nodes one by one
void saveQ(complex **Q, int ckpt_save) {
  int i, j, Ndat = 16 * DIMF, mynode, Nlines = 0, Qlen;
  char outfile[80];
  FILE *fp = NULL;

  // Total length of each column of Q
  Qlen = sites_on_node * Ndat * sizeof(complex);

  if (this_node == 0) {
    printf("Dumping columns %d--%d\n", ckpt_save + 1, volume * Ndat);
    sprintf(outfile, "%s.Q%d", startfile, ckpt_save);
    fp = fopen(outfile, "w");    // Open to write
    if (fp == NULL) {
      printf("saveQ: node%d can't open file %s\n", this_node, outfile);
      fflush(stdout);
      terminate(1);
    }

    // node0 data ready to go
    for (i = ckpt_save; i < volume * Ndat; i++) {
      for (j = 0; j < sites_on_node * Ndat; j++)
        if (Q[i][j].real != 0.0 || Q[i][j].real != 0.0) { // Loose condition
          fprintf(fp, "0\t%d\t%d\t%.16g\t%.16g\n",
                      i, j, Q[i][j].real, Q[i][j].imag);
          Nlines++;
      }
    }
    fflush(fp);
  }
  g_sync();   // Don't let other nodes race ahead

  // Overwrite node0 Q with elements of Q on mynode
  for (mynode = 1; mynode < numnodes(); mynode++) {
    for (i = ckpt_save; i < volume * Ndat; i++) {
      if (this_node == mynode)
        send_field((char *)Q[i], Qlen, 0);
      else if (this_node == 0)
        get_field((char *)Q[i], Qlen, mynode);

      g_sync();   // To be cautious
    }
    if (this_node == 0) {   // Now print, as above
      for (i = ckpt_save; i < volume * Ndat; i++) {
        for (j = 0; j < sites_on_node * Ndat; j++)
          if (Q[i][j].real != 0.0 || Q[i][j].real != 0.0) { // Loose condition
            fprintf(fp, "%d\t%d\t%d\t%.16g\t%.16g\n",
                        mynode, i, j, Q[i][j].real, Q[i][j].imag);
            Nlines++;
        }
      }
      fflush(fp);
    }
    g_sync();   // Don't let other nodes race ahead
  }

  if (this_node == 0) {
    // Add a last line containing a unique tag (-1)
    // and the number of elements that were written,
    // so we can check that the correct number were reloaded
    fprintf(fp, "-1\t%d\t-1\t-1\t-1\n", Nlines);
    fflush(fp);
    fclose(fp);
  }
  g_sync();   // Don't let other nodes race ahead
}
#endif
// -----------------------------------------------------------------
