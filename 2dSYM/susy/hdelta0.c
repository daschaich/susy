// -----------------------------------------------------------------
// Twist_Fermion matrix--vector operation, D^2 + fmass^2
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest = (D^2 + fmass^2).src
void hdelta0_field(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real fmass2 = fmass * fmass;
  Twist_Fermion *tmpTF;

  tmpTF = (Twist_Fermion *)malloc(sites_on_node * sizeof(Twist_Fermion));
  fermion_op(src, tmpTF, PLUS);
  fermion_op(tmpTF, dest, MINUS);
  if (fmass2 > 1e-6) {
    FORALLSITES(i, s)
      scalar_mult_add_TF(&(dest[i]), &(src[i]), fmass2, &(dest[i]));
  }
  free(tmpTF);
}
// -----------------------------------------------------------------
