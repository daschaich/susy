// -----------------------------------------------------------------
// Routines for susy gauge configuration input/output
// Loop over mu = 0 to NUMLINK
// Works for most machines
// Wrappers for I/O are in io_ansi.c

#include "generic_includes.h"
#include "../include/io_lat.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>

// This is dangerous -- it assumes off_t = long for this compilation
#ifndef HAVE_FSEEKO
#define fseeko fseek
#endif

#define EPS 1e-6

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* Checksums
   The dataset from which each checksum is computed is the full gauge
   configuration for lattice files and for propagator files, the
   propagator for a single source spin-color combination.  Data in these
   files appear as a series of 32-bit floating point numbers.  We treat
   the 32-bit values as unsigned integers v(i) where i = 0, ..., N-1 ranges
   over the values in the order of appearance on the file and N is the
   total count of values in the data set.  The checksum is obtained by a
   combination of bit-rotations and exclusive or operations.  It is
   designed to be commutative and associative, unlike the BSD sum
   operation, so that the data set can be read in parallel with checksum
   contributions computed for each portion read, and then combined
   afterwards.  The sum29 checksum does a left bit rotation through i mod
   29 bits and forms an exclusive or with the accumulated checksum.  The
   sum31 checksum does the same thing, but with i mod 31 bits.

   In writing the file the bit rotation is done on the number as
   represented on the architecture and consequently as written on the
   file.  In reading and checking file integrity on an architecure with a
   relatively byte-reversed representation, byte reversal of the data
   must be done before doing the bit rotation and the resulting checksum
   must be compared with the checksum recorded on the file after
   byte-reversal.
*/
// For checksums we want a 32 bit unsigned int, for which
// we have u_int32type defined in include/int32type.h which is
// included in include/io_lat.h

#define SUCCESS  0
#define FAILURE -1
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512

#define TOL 0.0000001   // Tolerance for floating point checks
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy a single precision fundamental matrix to generic precision
void f2d_mat(fmatrix *a, matrix *b) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      set_complex_equal(&(a->e[i][j]), &(b->e[i][j]));
  }
}

// Copy a generic precision fundamental matrix to single precision
void d2f_mat(matrix *a, fmatrix *b) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      set_complex_equal(&(a->e[i][j]), &(b->e[i][j]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Open a binary file for serial writing by node 0
// Return a file structure describing the opened file
gauge_file *w_serial_i(char *filename) {
  FILE *fp;
  gauge_file *gf;
  gauge_header *gh;

  // Set up gauge file and gauge header structures and load header values
  gf = setup_output_gauge_file();
  gh = gf->header;

  // Set number of nodes to zero to indicate coordinate natural ordering
  gh->order = NATURAL_ORDER;

  // Only node 0 opens the requested file
  if (this_node == 0) {
    fp = fopen(filename, "wb");
    if (fp == NULL) {
      printf("w_serial_i: node%d can't open file %s, error %d\n",
             this_node, filename, errno);
      fflush(stdout);
      terminate(1);
    }
    swrite_gauge_hdr(fp, gh);   // Write the header
  }

  // Assign values to file structure
  if (this_node == 0)   // Only node 0 knows about this file
    gf->fp = fp;
  else
    gf->fp = NULL;

  gf->filename = filename;
  gf->byterevflag = 0;    // Not used for writing
  gf->rank2rcv = NULL;    // Not used for writing
  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Flush lbuf to output, resetting buf_length
static void flush_lbuf_to_file(gauge_file *gf, fmatrix *lbuf,
                               int *buf_length) {

  FILE *fp = gf->fp;
  int stat;

  if (*buf_length <= 0)
    return;

  stat = (int)fwrite(lbuf, (NUMLINK + 2) * sizeof(fmatrix), *buf_length, fp);
  if (stat != *buf_length) {
    printf("w_serial: node%d gauge configuration write error %d file %s\n",
           this_node, errno, gf->filename);
    fflush(stdout);
    terminate(1);
  }
  *buf_length = 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate checksums
static void accum_cksums(gauge_file *gf, int *rank29, int *rank31,
                         u_int32type *buf, int n) {

  int k;
  u_int32type *val;

  for (k = 0, val = buf; k < n; k++, val++) {
    gf->check.sum29 ^= (*val)<<(*rank29) | (*val)>>(32-(*rank29));
    gf->check.sum31 ^= (*val)<<(*rank31) | (*val)>>(32-(*rank31));
    (*rank29)++;
    if (*rank29 >= 29)
      *rank29 = 0;
    (*rank31)++;
    if (*rank31 >= 31)
      *rank31 = 0;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Flush tbuf to lbuf and accumulate checksums without resetting tbuf_length
static void flush_tbuf_to_lbuf(gauge_file *gf, int *rank29, int *rank31,
                               fmatrix *lbuf, int *buf_length,
                               fmatrix *tbuf, int tbuf_length) {

  int nword;
  u_int32type *buf;

  if (tbuf_length > 0) {
    memcpy((void *)&lbuf[(NUMLINK + 2) * (*buf_length)],
           (void *)tbuf, (NUMLINK + 2) * tbuf_length * sizeof(fmatrix));

    nword = (NUMLINK + 2) * (int)sizeof(fmatrix)
                    / (int)sizeof(int32type) * tbuf_length;
    buf = (u_int32type *)&lbuf[(NUMLINK + 2) * (*buf_length)];
    accum_cksums(gf, rank29, rank31, buf, nword);

    *buf_length += tbuf_length;
  }
}

static void send_buf_to_node0(fmatrix *tbuf, int tbuf_length,
                              int currentnode) {

  if (this_node == currentnode) {
    send_field((char *)tbuf,
               (NUMLINK + 2) * tbuf_length * sizeof(fmatrix), 0);
  }
  else if (this_node == 0) {
    get_field((char *)tbuf,
              (NUMLINK + 2) * tbuf_length * sizeof(fmatrix), currentnode);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Only node 0 writes the gauge configuration to a binary file gf
void w_serial(gauge_file *gf) {
  register int i;
  int rank29, rank31, buf_length, tbuf_length, index;
  int j,x, t, currentnode, newnode;
  FILE *fp = NULL;
  gauge_header *gh = NULL;
  fmatrix *lbuf = NULL;
  fmatrix *tbuf = malloc(sizeof *tbuf * nx * (NUMLINK + 2));
  off_t offset;               // File stream pointer
  off_t coord_list_size;      // Size of coordinate list in bytes
  off_t head_size;            // Size of header plus coordinate list
  off_t checksum_offset = 0;  // Location of checksum
  off_t gauge_check_size;     // Size of checksum record

  // tbuf holds message buffer space for the x dimension
  // of the local hypercube (needs at most nx * (NUMLINK + 2) matrices)
  if (tbuf == NULL) {
    printf("w_serial: node%d can't malloc tbuf\n", this_node);
    terminate(1);
  }

  // Only allocate lbuf on node0
  if (this_node == 0) {
    lbuf = malloc(sizeof *lbuf * MAX_BUF_LENGTH * (NUMLINK + 2));
    if (lbuf == NULL) {
      printf("w_serial: node0 can't malloc lbuf\n");
      fflush(stdout);
      terminate(1);
    }

    fp = gf->fp;
    gh = gf->header;

    // No coordinate list written
    // Fields to be written in standard coordinate list order
    coord_list_size = 0;
    head_size = gh->header_bytes + coord_list_size;

    checksum_offset = head_size;

    gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);

    offset = head_size + gauge_check_size;

    if (fseeko(fp, offset, SEEK_SET) < 0) {
      printf("w_serial: node%d fseeko %lld failed error %d file %s\n",
             this_node, (long long)offset, errno, gf->filename);
      fflush(stdout);
      terminate(1);
    }
  }

  // Buffered algorithm for writing fields in serial (lexicographic) order
  // Initialize checksums
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;
  // Count 32-bit words mod 29 and mod 31 in order of appearance on file
  // Here only node 0 uses these values -- both start at 0
  i = sizeof(fmatrix) / sizeof(int32type) * sites_on_node * this_node;
  rank29 = ((NUMLINK + 2) * i) % 29;
  rank31 = ((NUMLINK + 2) * i) % 31;

  g_sync();
  currentnode = 0;  // The node delivering data
  buf_length = 0;
  tbuf_length = 0;
  for (t = 0; t < nt; t++) {
    for (x = 0; x < nx; x++) {
      newnode = node_number(x, t);  // The node providing the next site
      if (newnode != currentnode || x == 0) {
        // We are switching to a new node or have exhausted a line of nx
        // Sweep any data in the retiring node's tbuf to the node0 lbuf
        if (tbuf_length > 0) {
          if (currentnode != 0)
            send_buf_to_node0(tbuf, tbuf_length, currentnode);

          // node0 flushes tbuf, accumulates checksum
          // and writes lbuf if it is full
          if (this_node == 0) {
            flush_tbuf_to_lbuf(gf, &rank29, &rank31, lbuf, &buf_length,
                                                     tbuf, tbuf_length);
            if (buf_length > MAX_BUF_LENGTH - nx)
              flush_lbuf_to_file(gf, lbuf, &buf_length);
          }
          tbuf_length = 0;
        }

        // node0 sends a few bytes to newnode as a clear to send signal
        if (newnode != currentnode) {
          if (this_node == 0 && newnode != 0)
            send_field((char *)tbuf, NUMLINK, newnode);
          if (this_node == newnode && newnode != 0)
            get_field((char *)tbuf, NUMLINK, 0);
          currentnode = newnode;
        }
      }

      // The node with the data just appends to its tbuf
      if (this_node == currentnode) {
        i = node_index(x, t);
        index = (NUMLINK + 2) * tbuf_length;
        for (j = 0; j < NUMLINK; j++) {
          d2f_mat(&lattice[i].link[j], &tbuf[index]);
          index++;
        }
        d2f_mat(&lattice[i].phi, &tbuf[index]);
        index++;
        d2f_mat(&lattice[i].varphi, &tbuf[index]);
      }

      // tbuf_length is in units of (NUMLINK + 2) matrices
      if (this_node == currentnode || this_node == 0)
        tbuf_length++;
    }
  }

  // Purge any remaining data
  if (tbuf_length > 0) {
    if (currentnode != 0)
      send_buf_to_node0(tbuf, tbuf_length, currentnode);
  }

  if (this_node == 0) {
    flush_tbuf_to_lbuf(gf, &rank29, &rank31, lbuf, &buf_length,
                       tbuf, tbuf_length);
    flush_lbuf_to_file(gf, lbuf, &buf_length);
  }

  g_sync();
  free(tbuf);

  if (this_node == 0) {
    free(lbuf);
    printf("Saved gauge configuration serially to binary file %s\n",
           gf->filename);
    printf("Time stamp %s\n", gh->time_stamp);

    // Write checksum
    // Position file pointer
    if (fseeko(fp, checksum_offset, SEEK_SET) < 0) {
      printf("w_serial: node%d fseeko %lld failed error %d file %s\n",
             this_node, (long long)checksum_offset, errno, gf->filename);
      fflush(stdout);
      terminate(1);
    }
    write_checksum(gf);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Only node 0 reads the gauge configuration gf from a binary file
void r_serial(gauge_file *gf) {
  FILE *fp = gf->fp;
  gauge_header *gh = gf->header;
  char *filename = gf->filename;
  int byterevflag = gf->byterevflag;

  off_t offset = 0;           // File stream pointer
  off_t gauge_check_size;     // Size of gauge configuration checksum record
  off_t coord_list_size;      // Size of coordinate list in bytes
  off_t head_size;            // Size of header plus coordinate list
  off_t checksum_offset = 0;  // Where we put the checksum
  int rcv_rank, rcv_coords, destnode, stat, idest = 0;
  int j, k, x, t;
  int buf_length = 0, where_in_buf = 0;
  gauge_check test_gc;
  u_int32type *val;
  int rank29 = 0, rank31 = 0;
  fmatrix *lbuf = NULL;   // Only allocate on node0
  fmatrix tmat[NUMLINK + 2];

  if (this_node == 0) {
    // Compute offset for reading gauge configuration
    if (gh->magic_number == GAUGE_VERSION_NUMBER)
      gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);
    else
      gauge_check_size = 0;

    if (gf->header->order == NATURAL_ORDER)
      coord_list_size = 0;
    else
      coord_list_size = sizeof(int32type) * volume;
    checksum_offset = gf->header->header_bytes + coord_list_size;
    head_size = checksum_offset + gauge_check_size;

    // Allocate single-precision read buffer
    lbuf = malloc(sizeof *lbuf * MAX_BUF_LENGTH * (NUMLINK + 2));
    if (lbuf == NULL) {
      printf("r_serial: node%d can't malloc lbuf\n", this_node);
      fflush(stdout);
      terminate(1);
    }

    /* Position file for reading gauge configuration */
    offset = head_size;

    if (fseeko(fp, offset, SEEK_SET) < 0) {
      printf("r_serial: node0 fseeko %lld failed error %d file %s\n",
             (long long)offset, errno, filename);
      fflush(stdout);
      terminate(1);
    }
    buf_length = 0;
    where_in_buf = 0;
  }

  // All nodes initialize checksums
  // Will count 32-bit words mod 29 and mod 31 in order of appearance in file
  // All nodes see the same sequence because we read serially
  test_gc.sum29 = 0;
  test_gc.sum31 = 0;

  g_sync();

  // node0 reads and deals out the values
  for (rcv_rank = 0; rcv_rank < volume; rcv_rank++) {
    /* If file is in coordinate natural order, receiving coordinate
       is given by rank. Otherwise, it is found in the table */
    if (gf->header->order == NATURAL_ORDER)
      rcv_coords = rcv_rank;
    else
      rcv_coords = gf->rank2rcv[rcv_rank];

    x = rcv_coords % nx;
    rcv_coords /= nx;
    t = rcv_coords % nt;

    // The node that gets the next set of gauge links
    destnode = node_number(x, t);

    // node0 fills its buffer, if necessary
    if (this_node == 0) {
      if (where_in_buf == buf_length) {  /* get new buffer */
        /* new buffer length  = remaining sites, but never bigger
           than MAX_BUF_LENGTH */
        buf_length = volume - rcv_rank;
        if (buf_length > MAX_BUF_LENGTH)
          buf_length = MAX_BUF_LENGTH;

        // Now do read
        stat = (int)fread(lbuf, (NUMLINK + 2) * sizeof(fmatrix), buf_length, fp);
        if (stat != buf_length) {
          printf("r_serial: node%d gauge configuration read error %d file %s\n",
                 this_node, errno, filename);
          fflush(stdout);
          terminate(1);
        }
        where_in_buf = 0;  // Reset counter
      }  // End of the buffer read

      if (destnode == 0) {  // Just copy links
        idest = node_index(x, t);
        // Save (NUMLINK + 2) matrices in tmat for further processing
        memcpy(tmat, &lbuf[(NUMLINK + 2) * where_in_buf],
               (NUMLINK + 2) * sizeof(fmatrix));
      }
      else {                // Send to correct node
        send_field((char *)&lbuf[(NUMLINK + 2) * where_in_buf],
                   (NUMLINK + 2) * sizeof(fmatrix), destnode);
      }
      where_in_buf++;
    }

    // The node that contains this site reads the message
    else {  // All nodes other than node 0
      if (this_node == destnode) {
        idest = node_index(x, t);
        // Receive (NUMLINK + 2) matrices in temporary space for further processing
        get_field((char *)tmat, (NUMLINK + 2) * sizeof(fmatrix), 0);
      }
    }

    /* The receiving node does the byte reversal and then checksum,
       if needed.  At this point tmat contains the input matrices
       and idest points to the destination site structure. */
    if (this_node == destnode) {
      if (byterevflag == 1)
        byterevn((int32type *)tmat,
                 (NUMLINK + 2) * sizeof(fmatrix) / sizeof(int32type));
      // Accumulate checksums
      for (k = 0, val = (u_int32type *)tmat;
           k < (NUMLINK + 2) * (int)sizeof(fmatrix) / (int)sizeof(int32type);
           k++, val++) {
        test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32 - rank29);
        test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32 - rank31);
        rank29++;
        if (rank29 >= 29)
          rank29 = 0;
        rank31++;
        if (rank31 >= 31)
          rank31 = 0;
      }
      // Copy (NUMLINK + 2) matrices to generic-precision lattice[idest]
      for (j = 0; j < NUMLINK; j++)
        f2d_mat(&tmat[j], &lattice[idest].link[j]);
      f2d_mat(&tmat[NUMLINK], &lattice[idest].phi);
      f2d_mat(&tmat[NUMLINK + 1], &lattice[idest].varphi);

    }
    else {
      rank29 += (NUMLINK + 2) * sizeof(fmatrix) / sizeof(int32type);
      rank31 += (NUMLINK + 2) * sizeof(fmatrix) / sizeof(int32type);
      rank29 %= 29;
      rank31 %= 31;
    }
  }
  // Combine node checksum contributions with global exclusive or
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);

  if (this_node == 0) {
    // Read and verify checksum
    printf("Restored binary gauge and scalar configuration serially ");
    printf("from file %s\n", filename);
    if (gh->magic_number == GAUGE_VERSION_NUMBER) {
      printf("Time stamp %s\n", gh->time_stamp);
      if (fseeko(fp, checksum_offset, SEEK_SET) < 0) {
        printf("r_serial: node0 fseeko %lld failed error %d file %s\n",
               (long long)offset, errno, filename);
        fflush(stdout);
        terminate(1);
      }
      read_checksum(gf, &test_gc);
    }
    fflush(stdout);
    free(lbuf);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Top level routines
// Restore lattice file by reading serially from node 0
// Handles most lattice formats
gauge_file* restore_serial(char *filename) {
  gauge_file *gf;
  gf = r_serial_i(filename);
  if (gf->header->magic_number == LIME_MAGIC_NO) {
    r_serial_f(gf);
    // Close this reader and die with an error
    free(gf->header);
    free(gf);
    node0_printf("Looks like a SciDAC file -- unsupported\n");
    terminate(1);
  }
  else {
    r_serial(gf);
    r_serial_f(gf);
  }
  return gf;
}

// Save lattice in natural order by writing serially from node 0
gauge_file* save_serial(char *filename) {
  gauge_file *gf;

  gf = w_serial_i(filename);
  w_serial(gf);
  w_serial_f(gf);

  return gf;
}
// -----------------------------------------------------------------
