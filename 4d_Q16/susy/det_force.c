// -----------------------------------------------------------------
// Update the momenta with the determinant force
// Use tempmat and staple for temporary storage
#include "susy_includes.h"

double det_force(Real eps) {
  register int i, a, b, mu, nu;
  register site *s;
  double returnit = 0.0;
  msg_tag *tag[NUMLINK];

#ifdef TRUNCATED
  node0_printf("ERROR: Do not use non-zero kappa_u1 ");
  node0_printf("with truncated action... aborting\n");
  terminate(1);
#endif

  // Loop over directions, update momenta by accumulating
  //   ZWstar[b][a](x) + ZWstar[a][b](x - b) in tr_dest(x)
  // where ZWstar[a][b](x) = plaqdet[a][b](x) [plaqdet[a][b](x) - 1]^*
  FORALLDIR(a) {
    FORALLSITES(i, s)
      tr_dest[i] = cmplx(0.0, 0.0);

    // Start first gather (a = 0 and b = 1), labelled by b
    // Gather ZWstar[a][b] from x - b
    tag[1] = start_gather_field(ZWstar[0][1], sizeof(complex),
                                goffset[1] + 1, EVENANDODD, gen_pt[1]);

    // Main loop
    FORALLDIR(b) {
      if (a == b)
        continue;

      if (a < NUMLINK - 1 || b < NUMLINK - 2) {     // Start next gather
        if (b == NUMLINK - 1) {
          mu = a + 1;
          nu = 0;
        }
        else if (b == a - 1) {
          mu = a;
          nu = b + 2;
        }
        else {
          mu = a;
          nu = b + 1;
        }
        tag[nu] = start_gather_field(ZWstar[mu][nu], sizeof(complex),
                                     goffset[nu] + 1, EVENANDODD, gen_pt[nu]);
      }

      // Add ZWstar[b][a](x) to sum while gather runs
      FORALLSITES(i, s)
        CSUM(tr_dest[i], ZWstar[b][a][i]);

      // Add ZWstar[a][b](x - b) to sum
      wait_gather(tag[b]);
      FORALLSITES(i, s)
        CSUM(tr_dest[i], *((complex *)(gen_pt[b][i])));
      cleanup_gather(tag[b]);
    }

    // Now update momenta, overwriting f_U
    // Take adjoint and subtract to reproduce -Adj(f_U)
    // Factor of 2 obtained by summing over all b != a above,
    // compared to just b > a in action...
    FORALLSITES(i, s) {
      CMULREAL(tr_dest[i], kappa_u1, tr_dest[i]);
      c_scalar_mult_mat(&(Uinv[a][i]), &(tr_dest[i]), &(s->f_U[a]));
      scalar_mult_dif_adj_matrix(&(s->f_U[a]), eps, &(s->mom[a]));
    }
  }

  // Compute contribution to average gauge force
  // This is combined with the usual gauge force terms in update_h.c
  FORALLSITES(i, s) {
    FORALLDIR(a)
      returnit += realtrace(&(s->f_U[a]), &(s->f_U[a]));
  }
  g_doublesum(&returnit);

  return returnit;
}
// -----------------------------------------------------------------
