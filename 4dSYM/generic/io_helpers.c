// -----------------------------------------------------------------
// General purpose IO-related routines for susy code
// Loop over directions up to NUMLINK and do not reunitarize
#include "generic_includes.h"
#include "../include/io_lat.h"
#include "../susy/susy_includes.h"    // For plaquette
// -----------------------------------------------------------------



// -----------------------------------------------------------------
gauge_file *save_lattice(int flag, char *filename) {
  double dtime;
  gauge_file *gf = NULL;

  plaquette(&g_ssplaq, &g_stplaq);
  d_linktrsum(&linktrsum);
  nersc_checksum = nersc_cksum();

  dtime = -dclock();
  switch(flag) {
    case SAVE_SERIAL:
      gf = save_serial(filename);
      break;
    case FORGET:
      gf = NULL;
      break;
    default:
      node0_printf("\nsave_lattice: ERROR: unknown type for saving lattice\n");
      terminate(1);
  }
  dtime += dclock();
  if (flag != FORGET)
    node0_printf("Time to save = %e\n", dtime);
#if PRECISION == 1
  node0_printf("CHECK PLAQ: %e %e\n", g_ssplaq, g_stplaq);
  node0_printf("CHECK NERSC LINKTR: %e CKSUM: %x\n",
               linktrsum.real / (Real)NCOL, nersc_checksum);
#else
  // Double precision
  node0_printf("CHECK PLAQ: %.16e %.16e\n", g_ssplaq, g_stplaq);
  node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
               linktrsum.real / (Real)NCOL, nersc_checksum);
#endif
  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set linkf to unit matrices
void coldlat() {
  register int i, j, k, dir;
  register site *s;

  for (dir = 0; dir < NUMLINK; dir++) {
    FORALLSITES(i, s) {
      for (j = 0; j < NCOL; j++) {
        for (k = 0; k < NCOL; k++) {
          if (j != k)
            s->linkf[dir].e[j][k] = cmplx(0.0, 0.0);
          else
            s->linkf[dir].e[j][k] = cmplx(1.0, 0.0);
        }
      }
    }
  }
  node0_printf("unit gauge NUMLINK configuration loaded\n");
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set linkf to funny matrices for debugging
void funnylat() {
  register int i, j, k, dir;
  register site *s;

  FORALLSITES(i, s) {
    for (dir = XUP; dir <= TUP; dir++) {
      s->linkf[dir].e[0][0] = cmplx((double)dir, (double)dir);
      for (j = 1; j < NCOL; ++j) {
        for (k = 1; k < NCOL; ++k)
          s->linkf[dir].e[j][k] = cmplx(10.0 * j * s->x, 10.0 * k * s->z);
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Reload a lattice in binary format, set to unit gauge or keep current
// Do not reunitarize
gauge_file *reload_lattice(int flag, char *filename) {
  double dtime;
  gauge_file *gf = NULL;

  dtime = -dclock();
  switch(flag) {
    case CONTINUE:        // Do nothing
      gf = NULL;
      break;
    case FRESH:           // Cold lattice
      coldlat();
      gf = NULL;
      break;
    case RELOAD_SERIAL:   // Read binary lattice serially
      gf = restore_serial(filename);
      break;
    default:
      node0_printf("reload_lattice: Bad startflag %d\n", flag);
      terminate(1);
  }
  dtime += dclock();
  if (flag != FRESH && flag != CONTINUE)
    node0_printf("Time to reload gauge configuration = %e\n", dtime);

  plaquette(&g_ssplaq, &g_stplaq);
  d_linktrsum(&linktrsum);
  nersc_checksum = nersc_cksum();

#if PRECISION == 1
  node0_printf("CHECK PLAQ: %e %e\n", g_ssplaq, g_stplaq);
  node0_printf("CHECK NERSC LINKTR: %e CKSUM: %x\n",
               linktrsum.real / (Real)NCOL, nersc_checksum);
#else             // Double precision
  node0_printf("CHECK PLAQ: %.16e %.16e\n", g_ssplaq, g_stplaq);
  node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
               linktrsum.real / (Real)NCOL, nersc_checksum);
#endif
  fflush(stdout);
  dtime = -dclock();
  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Find out kind of starting lattice to use, and lattice name if necessary
// This routine is only called by node 0
int ask_starting_lattice(FILE *fp, int prompt, int *flag, char *filename) {
  char savebuf[256];
  int status;

  if (prompt!=0)
    printf("enter 'continue', 'fresh' or 'reload_serial'\n");
  status = fscanf(fp, "%s", savebuf);
  if (status == EOF) {
    printf("ask_starting_lattice: EOF on STDIN.\n");
    return 1;
  }
  if (status != 1) {
    printf("\nask_starting_lattice: ERROR IN INPUT: ");
    printf("can't read starting lattice option\n");
    return 1;
  }

  printf("%s", savebuf);
  if (strcmp("fresh", savebuf) == 0) {
    *flag = FRESH;
    printf("\n");
  }
  else if (strcmp("continue", savebuf) == 0) {
    *flag = CONTINUE;
    printf("\n");
  }
  else if (strcmp("reload_serial", savebuf) == 0)
    *flag = RELOAD_SERIAL;
  else {
    printf(" is not a valid starting lattice option. INPUT ERROR.\n");
    return 1;
  }

  // Read name of file and load it
  if (*flag != FRESH && *flag != CONTINUE) {
    if (prompt != 0)
      printf("enter name of file containing lattice\n");
    status = fscanf(fp, " %s", filename);
    if (status != 1) {
      printf("\nask_starting_lattice: ERROR IN INPUT: ");
      printf("error reading file name\n");
      return 1;
    }
    printf(" %s\n", filename);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Find out what do to with lattice at end, and lattice name if necessary
// This routine is only called by node 0
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename) {
  char savebuf[256];
  int status;

  if (prompt!=0)
    printf("'forget' lattice at end or 'save_serial'\n");
  status = fscanf(fp,"%s", savebuf);
  if (status != 1) {
    printf("\nask_ending_lattice: ERROR IN INPUT: error reading ending lattice command\n");
    return 1;
  }
  printf("%s", savebuf);
  if (strcmp("save_serial", savebuf) == 0)
    *flag = SAVE_SERIAL;
  else if (strcmp("forget", savebuf) == 0) {
    *flag = FORGET;
    printf("\n");
  }
  else {
    printf(" is not a save lattice command. INPUT ERROR\n");
    return 1;
  }

  if (*flag != FORGET) {
    if (prompt != 0)
      printf("enter filename\n");
    status = fscanf(fp, "%s", filename);
    if (status != 1) {
      printf("\nask_ending_lattice: ERROR IN INPUT: ");
      printf("error reading filename\n");
      return 1;
    }
    printf(" %s\n", filename);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Find out whether to gauge fix to Coulomb gauge
int ask_gauge_fix(FILE *fp, int prompt, int *flag) {
  int status = 0;
  char savebuf[256];

  if (prompt != 0)
    printf("enter 'no_gauge_fix', or 'coulomb_gauge_fix'\n");
  status = fscanf(fp, "%s", savebuf);
  if (status == EOF) {
    printf("ask_gauge_fix: EOF on STDIN\n");
    return 1;
  }
  if (status != 1) {
    printf("\nask_gauge_fix: ERROR IN INPUT: ");
    printf("can't read gauge fixing option\n");
    return 1;
  }

  printf("%s\n", savebuf);
  if (strcmp("coulomb_gauge_fix", savebuf) == 0)
    *flag = COULOMB_GAUGE_FIX;
  else if (strcmp("no_gauge_fix", savebuf) == 0)
    *flag = NO_GAUGE_FIX;
  else {
    printf("Error in input: invalid gauge fixing option\n");
    printf("Only no_gauge_fix and coulomb_gauge_fix supported\n");
    return 1;
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read and echo the next tag.  Echo any intervening comments
// Comments begin with # and apply to the rest of the line
// Verify that the input tag agrees with the expected tag
static int get_tag(FILE *fp, char *tag, char *myname) {
  static char checktag[80];
  char line[512];
  int s;

  while (1) {
    s = fscanf(fp, "%s", checktag);
    if (s == EOF) {
      printf("%s(%d): EOF on input\n", myname, this_node);
      return 1;
    }
    if (s == 0) {
      printf("%s(%d) Error reading %s\n", myname, this_node, tag);
      return 1;
    }
    if (strchr(checktag, '#') != NULL) {
      printf("%s", checktag);
      if (fgets(line, 512, fp) == NULL) {
        printf("%s(%d) EOF on input.\n", myname, this_node);
        return 1;
      }
      printf("%s", line);
    }
    else {
      if (strcmp(checktag, tag) != 0) {
        printf("\n%s: ERROR IN INPUT: expected %s but found %s\n",
               myname, tag, checktag);
        return 1;
      }
      printf("%s ", tag);
      return 0;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Check return value of scanf
static int check_read(int s, char *myname, char *tag) {
  if (s == EOF) {
    printf("\n%s: Expecting value for %s but found EOF\n", myname, tag);
    return 1;
  }
  else if (s == 0) {
    printf("\n%s: Format error reading value for %s\n", myname, tag);
    return 1;
  }
  else
    return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get a floating point number
// If prompt is non-zero, ask for the input value with tag
// If prompt is zero, require that tag precede the input value
// Return the value and exit on error
int get_f(FILE *fp, int prompt, char *tag, Real *value) {
  int s;
  char checkvalue[80];
  char myname[] = "get_f";

  if (prompt) {
    s = 0;
    while (s != 1) {
      printf("enter %s ", tag);
      fscanf(fp,"%s", checkvalue);
#if PRECISION == 1
      s = sscanf(checkvalue, "%e", value);
#else
      s = sscanf(checkvalue, "%le", value);
#endif
      if (s == EOF)
        return 1;
      if (s == 0)
        printf("Data format error.\n");
      else printf("%s %g\n", tag, *value);
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1)
      return 1;

#if PRECISION == 1
    s = fscanf(fp, "%e", value);
#else
    s = fscanf(fp, "%le", value);
#endif
    if (check_read(s, myname, tag) == 1)
      return 1;

    printf("%g\n", *value);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get an integer with same behavior as get_f
int get_i(FILE *fp, int prompt, char *tag, int *value) {
  int s;
  char checkvalue[80];
  char myname[] = "get_i";

  if (prompt) {
    s = 0;
    while (s != 1) {
      printf("enter %s ", tag);
      fscanf(fp, "%s", checkvalue);
      s=sscanf(checkvalue,"%d", value);
      if (s == EOF) return 1;
      if (s == 0) printf("Data format error\n");
      else printf("%s %d\n", tag, *value);
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1)
      return 1;

    s = fscanf(fp, "%d", value);
    if (check_read(s, myname, tag) == 1)
      return 1;
    printf("%d\n", *value);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read a single word as a string with same behavior as get_f
int get_s(FILE *fp, int prompt, char *tag, char *value) {
  int s;
  char myname[] = "get_s";

  if (prompt) {
    s = 0;
    while (s != 1) {
      printf("enter %s ", tag);
      s = fscanf(fp, "%s", value);
      if (s == EOF) return 1;
      if (s == 0) printf("Data format error\n");
      else printf("%s %s\n", tag, value);
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1) return 1;

    s = fscanf(fp, "%s", value);
    if (check_read(s, myname, tag) == 1) return 1;
    printf("%s\n", value);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read a vector of integers with same behavior as get_f
int get_vi(FILE* fp, int prompt, char *tag, int *value, int nvalues) {
  int s, i;
  char myname[] = "get_vi";

  if (prompt) {
    s = 0;
    printf("enter %s with %d values", tag, nvalues);
    for (i = 0; i < nvalues; i++) {
      while (s != 1) {
        printf("\n[%d] ", i);
        s = fscanf(fp, "%d", value + i);
        if (s == EOF) return 1;
        if (s == 0) printf("Data format error\n");
        printf("%s %d\n", tag, value[i]);
      }
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1) return 1;

    for (i = 0; i < nvalues - 1; i++) {
      s = fscanf(fp, "%d", value + i);
      if (check_read(s, myname, tag) == 1) return 1;
      printf("%d ", value[i]);
    }
    s = fscanf(fp, "%d", value + nvalues - 1);
    if (check_read(s, myname, tag) == 1) return 1;
    printf("%d\n", value[nvalues - 1]);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read a vector of reals with same behavior as get_f
int get_vf(FILE* fp, int prompt, char *tag, Real *value, int nvalues) {
  int s, i;
  char myname[] = "get_vf";

  if (prompt) {
    s = 0;
    printf("enter %s with %d values", tag, nvalues);
    for (i = 0; i < nvalues; i++) {
      while (s != 1) {
        printf("\n[%d] ", i);
#if PRECISION == 1
        s = scanf("%e", value + i);
#else
        s = scanf("%le", value + i);
#endif
        if (s == EOF)
          return 1;
        if (s == 0)
          printf("Data format error\n");

        printf("%s %g\n", tag, *(value + i));
      }
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1)
      return 1;

    for (i = 0; i < nvalues; i++) {
#if PRECISION == 1
      s = fscanf(fp, "%e", value + i);
#else
      s = fscanf(fp, "%le", value + i);
#endif
      if (check_read(s, myname, tag) == 1)
        return 1;
      printf("%g ", value[i]);
    }
    printf("\n");
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get the initial value of prompt:
// 0 for reading from file, 1 prompts for input from terminal
// Should be called only by node 0
// Return 0 if successful, 1 if failure
int get_prompt(FILE *fp, int *prompt) {
  char initial_prompt[512];
  int status;
  char myname[] = "get_prompt";

  *prompt = -1;
  printf("type 0 for no prompts or 1 for prompts\n");
  while (1) {
    status = fscanf(fp, "%s", initial_prompt);
    if (status != 1) {
      printf("\n%s: Can't read input\n", myname);
      terminate(1);
    }
    if (strchr(initial_prompt,'#') == NULL) break;
    // Provide for comment lines with # before "prompt"
    else {
      printf("%s", initial_prompt);
      if (fgets(initial_prompt, 512, fp) == NULL) {
        printf("%s(%d) EOF on input.\n", myname, this_node);
        return 1;
      }
      printf("%s", initial_prompt);
    }
  }
  if (strcmp(initial_prompt, "prompt") == 0)
    fscanf(fp, "%d", prompt);
  else if (strcmp(initial_prompt, "0") == 0)
    *prompt = 0;
  else if (strcmp(initial_prompt, "1") == 0)
    *prompt = 1;

  if (*prompt == 0 || *prompt == 1)
    return 0;
  else {
    printf("\n%s: ERROR IN INPUT: initial prompt\n", myname);
    return 1;
  }
}
// -----------------------------------------------------------------
