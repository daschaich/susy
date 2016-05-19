// -----------------------------------------------------------------
// Update the momenta with the determinant force
// Use tempmat and staple for temporary storage
#include "susy_includes.h"

double det_force(Real eps) {
  register int i, a, b;
  register site *s;
  double returnit = 0.0;
  msg_tag *mtag[8];

  // Loop over directions, update momenta
  FORALLDIR(a) {
    FORALLSITES(i, s)
      tr_dest[i] = cmplx(0.0, 0.0);

    // Loop over other directions,
    // accumulating ZWstar[b][a](x) + ZWstar[a][b](x - b) in tr_dest(x)
    // where ZWstar[a][b](x) = plaqdet[a][b](x) [plaqdet[a][b](x) - 1]^*
    FORALLDIR(b) {
      if (a == b)
        continue;

      mtag[0] = start_gather_field(ZWstar[a][b], sizeof(complex),
                                   goffset[b] + 1, EVENANDODD, gen_pt[0]);

      // Add ZWstar[b][a](x) to sum while gather runs
      FORALLSITES(i, s)
        CSUM(tr_dest[i], ZWstar[b][a][i]);

      // Add D[mu][nu](x - nu) to sum
      wait_gather(mtag[0]);
      FORALLSITES(i, s)
        CSUM(tr_dest[i], *((complex *)(gen_pt[0][i])));
      cleanup_gather(mtag[0]);
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
