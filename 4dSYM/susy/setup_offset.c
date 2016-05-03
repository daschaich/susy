// -----------------------------------------------------------------
// label[k] labels the offset of the kth connection between psibar and psi
// by an integer 1 to 4
// The separation of the paths is in New York metric
// It corresponds to our labeling lambda_1 = lambda(1, 0, 0, 0)
// offset[k] is the actual offset of the kth connection
// This program only generates the positive offsets, 40 of the 80 total
// The others can be found by using adjoints of these paths
// (e.g., (-1, -1, -1, -1) is the adjoint of (1, 1, 1, 1))

#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void cubic_neighbor(int x, int y, int z, int t, int *arg, int forw_back,
                    int *xpt, int *ypt, int *zpt, int *tpt) {

  if (forw_back == FORWARDS) {
    *xpt = (x + nx + arg[0]) % nx;
    *ypt = (y + ny + arg[1]) % ny;
    *zpt = (z + nz + arg[2]) % nz;
    *tpt = (t + nt + arg[3]) % nt;
  }
  else {
    *xpt = (x + nx - arg[0]) % nx;
    *ypt = (y + ny - arg[1]) % ny;
    *zpt = (z + nz - arg[2]) % nz;
    *tpt = (t + nt - arg[3]) % nt;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int check_list(int my_offset[4], int *minus) {
  int n, index = -1;
  *minus = 1;

  for (n = 0; n < q_off_max && index == -1; n++) {
    if (my_offset[0] == q_offset[n][0] && my_offset[1] == q_offset[n][1] &&
        my_offset[2] == q_offset[n][2] && my_offset[3] == q_offset[n][3])
      index = n;
  }

  // Check negative if needed
  if (index == -1) {
    *minus = -1;
    for (n = 0; n < q_off_max && index == -1; n++) {
      if (my_offset[0] == -q_offset[n][0] && my_offset[1] == -q_offset[n][1]
       && my_offset[2] == -q_offset[n][2] && my_offset[3] == -q_offset[n][3])
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
  register int i, dir1, dir2, dir3, tocheck;
  register site *s;

  // Single-offset terms
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    FORALLSITES(i, s) {
      s->bc1[dir1] = 1.0;
      s->bc1[OPP_LDIR(dir1)] = 1.0;

      if (s->t + offset[dir1][TUP] < 0)
        s->bc1[dir1] = PBC;
      else if (s->t + offset[dir1][TUP] > nt - 1)
        s->bc1[dir1] = PBC;
      if (s->t - offset[dir1][TUP] < 0)
        s->bc1[OPP_LDIR(dir1)] = PBC;
      else if (s->t - offset[dir1][TUP] > nt - 1)
        s->bc1[OPP_LDIR(dir1)] = PBC;
    }
  }

  // BC1 test
//  FORALLSITES(i, s) {
//    if (s->x == 1 && s->y == 1 && s->z == 2) {
//      printf("%d %d %d %d:", s->x, s->y, s->z, s->t);
//      for (dir1 = 0; dir1 < 2 * NUMLINK; dir1++)
//        printf("  %4.2f", s->bc1[dir1]);
//      printf("\n");
//    }
//  }

  // Double-offset terms -- don't need mixed -+ and +-
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    for (dir2 = 0; dir2 < NUMLINK; dir2++) {
      FORALLSITES(i, s) {
        s->bc2[dir1][dir2] = 1.0;
        s->bc2[OPP_LDIR(dir1)][OPP_LDIR(dir2)] = 1.0;

        tocheck = s->t + offset[dir1][TUP] + offset[dir2][TUP];
        if (tocheck < 0)
          s->bc2[dir1][dir2] = PBC;
        else if (tocheck > nt - 1)
          s->bc2[dir1][dir2] = PBC;

        tocheck = s->t - offset[dir1][TUP] - offset[dir2][TUP];
        if (tocheck < 0)
          s->bc2[OPP_LDIR(dir1)][OPP_LDIR(dir2)] = PBC;
        else if (tocheck > nt - 1)
          s->bc2[OPP_LDIR(dir1)][OPP_LDIR(dir2)] = PBC;
      }
    }
  }

  // Triple-offset terms -- don't need mixed -++, +-+, ++-, --+, -+- or +--
  for (dir1 = 0; dir1 < NUMLINK; dir1++) {
    for (dir2 = 0; dir2 < NUMLINK; dir2++) {
      for (dir3 = 0; dir3 < NUMLINK; dir3++) {
        FORALLSITES(i, s) {
          s->bc3[dir1][dir2][dir3] = 1.0;
          s->bc3[OPP_LDIR(dir1)][OPP_LDIR(dir2)][OPP_LDIR(dir3)] = 1.0;

          tocheck = s->t + offset[dir1][TUP] + offset[dir2][TUP]
                         + offset[dir3][TUP];
          if (tocheck < 0)
            s->bc3[dir1][dir2][dir3] = PBC;
          else if (tocheck > nt - 1)
            s->bc3[dir1][dir2][dir3] = PBC;

          tocheck = s->t - offset[dir1][TUP] - offset[dir2][TUP]
                         - offset[dir3][TUP];
          if (tocheck < 0)
            s->bc3[OPP_LDIR(dir1)][OPP_LDIR(dir2)][OPP_LDIR(dir3)] = PBC;
          else if (tocheck > nt - 1)
            s->bc3[OPP_LDIR(dir1)][OPP_LDIR(dir2)][OPP_LDIR(dir3)] = PBC;
        }
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_offset() {
  int i, k;

  // Construct the link paths: one in each direction plus back-diagonal
  for (i = 0; i < NUMLINK - 1; i++) {
    for (k = 0; k < NDIMS; k++)
      offset[i][k] = 0;

    offset[i][i] = 1;
  }
  for (k = 0; k < NDIMS; k++)
    offset[DIR_5][k] = -1;

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
  for (i = 0; i < NUMLINK; i++) {
    goffset[i] = make_gather(cubic_neighbor, offset[i],
                             WANT_INVERSE, NO_EVEN_ODD, SCRAMBLE_PARITY);

#ifdef DEBUG_CHECK
    int dir;
    node0_printf("  %d ahead:", i);
    for (dir = XUP; dir <= TUP; dir++)
      node0_printf(" %d", offset[i][dir]);

    node0_printf(" (offset %d)\n", goffset[i]);
#endif
  }

  // Make boundary condition tables
  if (PBC >= 0)
    node0_printf("Periodic temporal boundary conditions\n");
  if (PBC < 0)
    node0_printf("Antiperiodic temporal boundary conditions\n");
  setup_bc();
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void setup_qclosed_offset() {
  int flag, i, ii, mu;
  int a, b, c, d, e;
  int DbmP, DbpP, F1Q, F2Q;

  // Info for general gathers
  int n[4], k, d1[4], d2[4], d3[4];
  int ipath, d1_list, d2_list, minus;
  static int mu_vec[NUMLINK][NDIMS] = {{ 1,  0,  0,  0},
                                       { 0,  1,  0,  0},
                                       { 0,  0,  1,  0},
                                       { 0,  0,  0,  1},
                                       {-1, -1, -1, -1}};


  // Construct the first list of offsets
  // We only need those with 2 and 3 nonzero offsets
  q_off_max = 0;
  for (n[0] = 1; n[0] >= -1; n[0]--) {
    for (n[1] = 1; n[1] >= -1; n[1]--) {
      for (n[2] = 1; n[2] >= -1; n[2]--) {
        for (n[3] = 1; n[3] >= -1; n[3]--) {
          k = abs(n[0]) + abs(n[1]) + abs(n[2]) + abs(n[3]);
          if (k == 2 || k == 3) {
            // Run the list of previous entries
            // to see if the new one is already here, reversed
            flag = 0;
            for (ii = 0; ii < q_off_max; ii++) {
              if (n[0] == -q_offset[ii][0]
               && n[1] == -q_offset[ii][1]
               && n[2] == -q_offset[ii][2]
               && n[3] == -q_offset[ii][3])
                flag = 1;
            }
            if (flag == 0) {
              for (mu = 0; mu < NDIMS; mu++)
                q_offset[q_off_max][mu] = n[mu];

              q_off_max++;
            }
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

  // Set up gathers for DbminusPtoP terms
  DbmP = 0;
  for (d = 0; d < NUMLINK; d++) {
    for (e = d + 1; e < NUMLINK; e++) {
      for (c = 0; c < NUMLINK; c++) {
        if (c == d || c == e)
          continue;

        for (a = 0; a < NUMLINK; a++) {
          if (a == c || a == d || a == e)
            continue;

          for (b = a + 1; b < NUMLINK; b++) {
            if (b == c || b == d || b == e)
              continue;

            // General gathers d1 = -e_a - e_b - e_c
            //                 d2 = -e_a - e_b
            for (i = 0; i < NDIMS; i++) {
              d2[i] = -mu_vec[a][i] - mu_vec[b][i];
              d1[i] = d2[i] - mu_vec[c][i];
            }

            d1_list = check_list(d1, &minus);
            if (minus > 0)
              DbmP_d1[DbmP] = gq_offset[d1_list];
            if (minus < 0)
              DbmP_d1[DbmP] = gq_offset[d1_list] + 1;

            d2_list = check_list(d2, &minus);
            if (minus > 0)
              DbmP_d2[DbmP] = gq_offset[d2_list];
            if (minus < 0)
              DbmP_d2[DbmP] = gq_offset[d2_list] + 1;

#ifdef DEBUG_CHECK
            node0_printf("DbminusPtoP general gather d1 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d1[0], d1[1], d1[2], d1[3],
                         d1_list, DbmP_d1[DbmP]);
            node0_printf("DbminusPtoP general gather d2 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d2[0], d2[1], d2[2], d2[3],
                         d2_list, DbmP_d2[DbmP]);
#endif
            DbmP++;
          }
        }
      }
    }
  }
#ifdef DEBUG_CHECK
  node0_printf("Total of %d DbmP terms\n", DbmP);
#endif

  // Set up gathers for DbplusPtoP terms
  DbpP = 0;
  for (a = 0; a < NUMLINK; a++) {
    for (b = a + 1; b < NUMLINK; b++) {
      for (c = 0;c<NUMLINK;c++) {
        if (c == a || c == b)
          continue;
        for (d = 0; d < NUMLINK; d++) {
          if (d == c || d == a || d == b)
            continue;
          for (e = d + 1; e < NUMLINK; e++) {
            if ((e==c)||(e==a)||(e==b))
              continue;
            // General gathers d1 = e_a + e_b
            //                 d2 = e_a + e_b + e_c
            for (i = 0; i < NDIMS; i++) {
              d1[i] = mu_vec[a][i] + mu_vec[b][i];
              d2[i] = d1[i] + mu_vec[c][i];
            }

            d1_list = check_list(d1, &minus);
            if (minus > 0)
              DbpP_d1[DbpP] = gq_offset[d1_list];
            if (minus < 0)
              DbpP_d1[DbpP] = gq_offset[d1_list] + 1;

            d2_list = check_list(d2, &minus);
            if (minus > 0)
              DbpP_d2[DbpP] = gq_offset[d2_list];
            if (minus < 0)
              DbpP_d2[DbpP] = gq_offset[d2_list] + 1;

#ifdef DEBUG_CHECK
            node0_printf("DbplusPtoP general gather d1 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d1[0], d1[1], d1[2], d1[3],
                         d1_list, DbpP_d1[DbpP]);
            node0_printf("DbplusPtoP general gather d2 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d2[0], d2[1], d2[2], d2[3],
                         d2_list, DbpP_d2[DbpP]);
#endif
            DbpP++;
          }
        }
      }
    }
  }
#ifdef DEBUG_CHECK
  node0_printf("Total of %d DbpP terms\n", DbpP);
#endif

  // Set up gathers for first Q-closed force piece
  F1Q = 0;
  for (c = 0; c < NUMLINK; c++) {
    for (d = 0; d < NUMLINK; d++) {
      if (d == c)
        continue;
      for (e = d + 1; e < NUMLINK; e++) {
        if (e == c)
          continue;
        for (a = 0; a < NUMLINK; a++) {
          if (a == d || a == e || a == c)
            continue;
          for (b = a + 1; b < NUMLINK; b++) {
            if (b == d || b == e || b == c)
              continue;
            // General gathers d3 = e_a + e_b
            //                 d2 = e_a + e_b + e_c
            // Also gather d1 = -e_a - e_b
            for (i = 0; i < NDIMS; i++) {
              d3[i] = mu_vec[a][i] + mu_vec[b][i];
              d2[i] = d3[i] + mu_vec[c][i];
              d1[i] = -d3[i];
            }

            d1_list = check_list(d1, &minus);
            if (minus > 0)
              F1Q_d1[F1Q] = gq_offset[d1_list];
            if (minus < 0)
              F1Q_d1[F1Q] = gq_offset[d1_list] + 1;

            d2_list = check_list(d2, &minus);
            if (minus > 0)
              F1Q_d2[F1Q] = gq_offset[d2_list];
            if (minus < 0)
              F1Q_d2[F1Q] = gq_offset[d2_list] + 1;

#ifdef DEBUG_CHECK
            node0_printf("Qclosed1 general gather d1 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d1[0], d1[1], d1[2], d1[3],
                         d1_list, F1Q_d1[F1Q]);
            node0_printf("Qclosed1 general gather d2 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d2[0], d2[1], d2[2], d2[3],
                         d2_list, F1Q_d2[F1Q]);
#endif
            F1Q++;
          }
        }
      }
    }
  }
#ifdef DEBUG_CHECK
  node0_printf("Total of %d F1Q terms\n", F1Q);
#endif

  // Set up gathers for second Q-closed force piece
  F2Q = 0;
  for (c = 0; c < NUMLINK; c++) {
    for (d = 0; d < NUMLINK; d++) {
      if (d == c)
        continue;
      for (e = d + 1; e < NUMLINK; e++) {
        if (e == c)
          continue;
        for (a = 0; a < NUMLINK; a++) {
          if (a == d || a == e || a == c)
            continue;
          for (b = a + 1; b < NUMLINK; b++) {
            if (b == d || b == e || b == c)
              continue;

            // General gathers d3 = e_a + e_b
            //                 d2 = e_a + e_b + e_c
            // Also gather d1 = -e_a - e_b
            for (i = 0; i < NDIMS; i++) {
              d3[i] = mu_vec[a][i] + mu_vec[b][i];
              d2[i] = d3[i] + mu_vec[c][i];
              d1[i] = -d3[i];
            }

            d1_list = check_list(d1,&minus);
            if (minus > 0)
              F2Q_d1[F2Q] = gq_offset[d1_list];
            if (minus < 0)
              F2Q_d1[F2Q] = gq_offset[d1_list] + 1;

            d2_list = check_list(d2,&minus);
            if (minus > 0)
              F2Q_d2[F2Q] = gq_offset[d2_list];
            if (minus < 0)
              F2Q_d2[F2Q] = gq_offset[d2_list] + 1;

#ifdef DEBUG_CHECK
            node0_printf("Qclosed2 general gather d1 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d1[0], d1[1], d1[2], d1[3],
                         d1_list, F2Q_d1[F2Q]);
            node0_printf("Qclosed2 general gather d2 ");
            node0_printf("(%d %d %d %d), %d, %d\n",
                         d2[0], d2[1], d2[2], d2[3],
                         d2_list, F2Q_d2[F2Q]);
#endif
            F2Q++;
          }
        }
      }
    }
  }
#ifdef DEBUG_CHECK
  node0_printf("Total of %d F2Q terms\n", F2Q);
#endif
}
// -----------------------------------------------------------------
