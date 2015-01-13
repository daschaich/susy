// -----------------------------------------------------------------
// Twist_Fermion matrix--vector operation, D^2 + fmass^2
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// dest = (D^2 + fmass^2).src
// Use tempTF for temporary storage
void hdelta0_field(Twist_Fermion *src, Twist_Fermion *dest) {
  register int i;
  register site *s;
  Real fmass2 = fmass * fmass;

  fermion_op(src, tempTF, PLUS);
  fermion_op(tempTF, dest, MINUS);
  if (fmass2 > IMAG_TOL) {
    FORALLSITES(i, s)
      scalar_mult_add_TF(&(dest[i]), &(src[i]), fmass2, &(dest[i]));
  }
}
// -----------------------------------------------------------------
