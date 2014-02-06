// -----------------------------------------------------------------
// Wrapper for creating link field 'link' in irrep of fermions
// from field 'linkf' in fundamental rep
// Calls single-site routine to perform the actual translation
// Keep single-site routine in separate file(s) for various reps
#include "susy_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void fermion_rep() {
  site *s;
  int mu, i;

  for (mu = 0; mu < NUMLINK; mu++){
    FORALLSITES(i, s) {
      make_fermion_rep_matrix(&(s->linkf[mu]), &(s->link[mu]));
//      node0_printf("Site %d mu %d\n", i, mu);
//      dumpmat(&(s->link[mu]));
    }
  }
}
// -----------------------------------------------------------------
