// -----------------------------------------------------------------
// label[k] labels the offset of the kth connection between psibar and psi
// by an integer 1 to 3
// The separation of the paths is in New York metric
// It corresponds to our labeling lambda_1 = lambda(1, 0, 0)
// offset[k] is the actual offset of the kth connection
// This program only generates the positive offsets, 40 of the 80 total
// The others can be found by using adjoints of these paths
// (e.g., (-1, -1, -1) is the adjoint of (1, 1, 1))

#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void cubic_neighbor(int x, int y, int t, int *arg, int forw_back,
                    int *xpt, int *ypt, int *tpt) {

  if (forw_back == FORWARDS) {
    *xpt = (x + nx + arg[0]) % nx;
    *ypt = (y + ny + arg[1]) % ny;
    *tpt = (t + nt + arg[2]) % nt;
  }
  else {
    *xpt = (x + nx - arg[0]) % nx;
    *ypt = (y + ny - arg[1]) % ny;
    *tpt = (t + nt - arg[2]) % nt;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int check_list(int my_offset[3], int *minus) {
  int n, index = -1;
  *minus = 1;

  for (n = 0; n < q_off_max && index == -1; n++) {
    if (my_offset[0] == q_offset[n][0] && my_offset[1] == q_offset[n][1] &&
        my_offset[2] == q_offset[n][2])
      index = n;
  }

  // Check negative if needed
  if (index == -1) {
    *minus = -1;
    for (n = 0; n < q_off_max && index == -1; n++) {
      if (my_offset[0] == -q_offset[n][0] && my_offset[1] == -q_offset[n][1]
       && my_offset[2] == -q_offset[n][2])
        index = n;
    }
  }
  if (index < 0)
    node0_printf("Couldn't find entry in list\n");

  return index;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_bc() {
  register int i, dir, dir2, dir3, tocheck;
  register site *s;

  // Single-offset terms
  FORALLDIR(dir) {
    FORALLSITES(i, s) {
      s->bc[dir] = 1.0;
      s->bc[OPP_LDIR(dir)] = 1.0;

      if (s->t + offset[dir][TUP] < 0)
        s->bc[dir] = PBC;
      else if (s->t + offset[dir][TUP] > nt - 1)
        s->bc[dir] = PBC;
      if (s->t - offset[dir][TUP] < 0)
        s->bc[OPP_LDIR(dir)] = PBC;
      else if (s->t - offset[dir][TUP] > nt - 1)
        s->bc[OPP_LDIR(dir)] = PBC;
    }
  }

  // BC1 test
//  FORALLSITES(i, s) {
//    if (s->x == 1 && s->y == 1 && s->z == 2) {
//      printf("%d %d %d %d:", s->x, s->y, s->z, s->t);
//      for (dir = 0; dir < 2 * NUMLINK; dir++)
//        printf("  %4.2g", s->bc[dir]);
//      printf("\n");
//    }
//  }

  // Double-offset terms -- don't need mixed -+ and +-
  FORALLDIR(dir) {
    FORALLDIR(dir2) {
      FORALLSITES(i, s) {
        s->bc2[dir][dir2] = 1.0;
        s->bc2[OPP_LDIR(dir)][OPP_LDIR(dir2)] = 1.0;

        tocheck = s->t + offset[dir][TUP] + offset[dir2][TUP];
        if (tocheck < 0)
          s->bc2[dir][dir2] = PBC;
        else if (tocheck > nt - 1)
          s->bc2[dir][dir2] = PBC;

        tocheck = s->t - offset[dir][TUP] - offset[dir2][TUP];
        if (tocheck < 0)
          s->bc2[OPP_LDIR(dir)][OPP_LDIR(dir2)] = PBC;
        else if (tocheck > nt - 1)
          s->bc2[OPP_LDIR(dir)][OPP_LDIR(dir2)] = PBC;
      }
    }
  }

  // Triple-offset terms -- don't need mixed -++, +-+, ++-, --+, -+- or +--
  FORALLDIR(dir) {
    FORALLDIR(dir2) {
      FORALLDIR(dir3) {
        FORALLSITES(i, s) {
          s->bc3[dir][dir2][dir3] = 1.0;
          s->bc3[OPP_LDIR(dir)][OPP_LDIR(dir2)][OPP_LDIR(dir3)] = 1.0;

          tocheck = s->t + offset[dir][TUP] + offset[dir2][TUP]
                         + offset[dir3][TUP];
          if (tocheck < 0)
            s->bc3[dir][dir2][dir3] = PBC;
          else if (tocheck > nt - 1)
            s->bc3[dir][dir2][dir3] = PBC;

          tocheck = s->t - offset[dir][TUP] - offset[dir2][TUP]
                         - offset[dir3][TUP];
          if (tocheck < 0)
            s->bc3[OPP_LDIR(dir)][OPP_LDIR(dir2)][OPP_LDIR(dir3)] = PBC;
          else if (tocheck > nt - 1)
            s->bc3[OPP_LDIR(dir)][OPP_LDIR(dir2)][OPP_LDIR(dir3)] = PBC;
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_offset() {
  int i, k;

  // Construct the link paths: one in each direction
  FORALLUPDIR(i) {
    for (k = 0; k < NDIMS; k++)
      offset[i][k] = 0;

    offset[i][i] = 1;
  }
  // don't need back diagonal in 3d

#ifdef DEBUG_CHECK
  node0_printf("There are %d distinct paths:\n", NUMLINK);
#endif
  // goffset holds indices of gather_array in ../generic/com_mpi.c
  // The first eight elements of gather_array are
  //   XUP, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN
  // in that order!
  // In order to use XDOWN = XUP + 1, etc., we make the next ten elements
  //   XUP, XDOWN, YUP, YDOWN, ZUP, ZDOWN, TUP, TDOWN, DIR_5, -DIR_5
  // Then goffset[0]=8, goffset[1]=10, ..., goffset[4]=16
  // But we can't use these in EVEN or ODD gathers!
  FORALLDIR(i) {
    goffset[i] = make_gather(cubic_neighbor, offset[i],
                             WANT_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY);

#ifdef DEBUG_CHECK
    int dir;
    node0_printf("  %d ahead:", i);
    FORALLUPDIR(dir)
      node0_printf(" %d", offset[i][dir]);

    node0_printf(" (offset %d)\n", goffset[i]);
#endif
  }

  // Make boundary condition tables
  if (PBC == 1) {               // Braces suppress compiler warning
    node0_printf("Periodic temporal boundary conditions\n");
  }
  else if (PBC == -1) {         // Braces suppress compiler warning
    node0_printf("Antiperiodic temporal boundary conditions\n");
  }
  else {
    node0_printf("Error: Periodic or antiperiodic BCs expected\n");
    terminate(1);
  }
  setup_bc();
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_qclosed_offset() {

  int flag, i, ii, mu;
  int a, b, c;
  int DbVP, FQ;

  // Info for general gathers
  int n[3], k, d1[3];
  int ipath, d1_list, minus;
  static int mu_vec[NUMLINK][NDIMS] = {{ 1,  0,  0},
                                       { 0,  1,  0},
                                       { 0,  0,  1}};


  // Construct the first list of offsets
  // We only need those with 2 and 3 nonzero offsets
  q_off_max = 0;
  for (n[0] = 1; n[0] >= -1; n[0]--) {
    for (n[1] = 1; n[1] >= -1; n[1]--) {
      for (n[2] = 1; n[2] >= -1; n[2]--) {
        k = abs(n[0]) + abs(n[1]) + abs(n[2]);
        if (k == 2 || k == 3) {
          // Run the list of previous entries
          // to see if the new one is already here, reversed
          flag = 0;
          for (ii = 0; ii < q_off_max; ii++) {
            if (n[0] == -q_offset[ii][0]
             && n[1] == -q_offset[ii][1]
             && n[2] == -q_offset[ii][2])
              flag = 1;
          }
          if (flag == 0) {
            FORALLUPDIR(mu)
              q_offset[q_off_max][mu] = n[mu];

            q_off_max++;
          }
        }
      }
    }
  }


#ifdef DEBUG_CHECK
  node0_printf("There are %d distinct Q-closed candidate paths:\n", q_off_max);
#endif
  // Now make the gather tables
  for (ipath = 0; ipath < q_off_max; ipath++) {
    gq_offset[ipath] = make_gather(cubic_neighbor, q_offset[ipath],
                                   WANT_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY);

#ifdef DEBUG_CHECK
    node0_printf("  %d ahead:", ipath);
    for (mu = 0; mu < NDIMS; mu++)
      node0_printf(" %d", q_offset[ipath][mu]);

    node0_printf(" (offset %d)\n", gq_offset[ipath]);
#endif
  }

  // Set up gathers for Db(minus/plus)VtoP terms
  DbVP = 0;
  FORALLDIR(a) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLDIR(c) {
        if (a == c || b == c)
          continue;
    // General gathers d1 = e_a + e_b
        FORALLUPDIR(i) {
          d1[i] = mu_vec[a][i] + mu_vec[b][i];
        }

        d1_list = check_list(d1, &minus);
        if (minus > 0)
          DbVP_d1[DbVP] = gq_offset[d1_list];
        if (minus < 0)
          DbVP_d1[DbVP] = gq_offset[d1_list] + 1;
#ifdef DEBUG_CHECK
        node0_printf("DbVtoP general gather d1 ");
        node0_printf("(%d %d %d), %d, %d\n",
                     d1[0], d1[1], d1[2],
                     d1_list, DbVP_d1[DbVP]);
#endif
      DbVP++;
      }
    }
  }
#ifdef DEBUG_CHECK
  node0_printf("Total of %d DbVP terms\n", DbVP);
#endif

  //Set up gathers for the Q-closed terms
  //Both FQ and F2Q use FQ_d1
  FQ = 0;
  FORALLDIR(a) {
    for (b = a + 1; b < NUMLINK; b++) {
      FORALLDIR(c) {
        if (a == c || b == c)
          continue;
    // General gathers d1 = - e_a - e_b
        FORALLUPDIR(i) {
          d1[i] = - mu_vec[a][i] - mu_vec[b][i];
        }

        d1_list = check_list(d1, &minus);
        if (minus > 0)
          FQ_d1[FQ] = gq_offset[d1_list];
        if (minus < 0)
          FQ_d1[FQ] = gq_offset[d1_list] + 1;
#ifdef DEBUG_CHECK
        node0_printf("FQ general gather d1 ");
        node0_printf("(%d %d %d), %d, %d\n",
                     d1[0], d1[1], d1[2],
                     d1_list, FQ_d1[FQ]);
#endif
      FQ++;
      }
    }
  }
#ifdef DEBUG_CHECK
  node0_printf("Total of %d FQ terms\n", FQ);
#endif
}
// -----------------------------------------------------------------
