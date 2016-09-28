// -----------------------------------------------------------------
// Z2 distributed random number
// This requires a random number generator, myrand(),
#include "../include/config.h"
#include "../include/susy.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// prn_pt is a pointer passed to myrand()
Real z2_rand_no(double_prn *prn_pt) {
   if(myrand(prn_pt)<0.5){
	   return (-1.0);
   }
   return (1.0);
}
// -----------------------------------------------------------------
