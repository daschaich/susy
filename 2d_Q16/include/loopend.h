// -----------------------------------------------------------------
// Redefinition of FORSOMEPARITY for slightly better looping
// macros.h must come first
#ifndef _LOOPEND_H
#define _LOOPEND_H

#include "../include/macros.h"
#define LOOPEND   // We take this as the rule for now
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#ifdef LOOPEND
#undef FORSOMEPARITY
#define FORSOMEPARITY(i,s,choice) \
{ register int loopend;  \
loopend= (choice)==EVEN ? even_sites_on_node : sites_on_node ; \
for( i=((choice)==ODD ? even_sites_on_node : 0 ), s= &(lattice[i]); \
i<loopend; i++,s++)
#define END_LOOP }
#else
#define END_LOOP  // Defined to be nothing
#endif

#endif
// -----------------------------------------------------------------
