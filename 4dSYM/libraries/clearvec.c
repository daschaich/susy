// -----------------------------------------------------------------
// Clear an irrep vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec(su3_vector *v) {
  v->c[0].real = 0.0;     v->c[0].imag = 0.0;
  v->c[1].real = 0.0;     v->c[1].imag = 0.0;
  v->c[2].real = 0.0;     v->c[2].imag = 0.0;
  v->c[3].real = 0.0;     v->c[3].imag = 0.0;
#if (DIMF > 4)
  v->c[4].real = 0.0;     v->c[4].imag = 0.0;
  v->c[5].real = 0.0;     v->c[5].imag = 0.0;
  v->c[6].real = 0.0;     v->c[6].imag = 0.0;
  v->c[7].real = 0.0;     v->c[7].imag = 0.0;
  v->c[8].real = 0.0;     v->c[8].imag = 0.0;
#if (DIMF > 9)
  v->c[9].real = 0.0;     v->c[9].imag = 0.0;
  v->c[10].real = 0.0;    v->c[10].imag = 0.0;
  v->c[11].real = 0.0;    v->c[11].imag = 0.0;
  v->c[12].real = 0.0;    v->c[12].imag = 0.0;
  v->c[13].real = 0.0;    v->c[13].imag = 0.0;
  v->c[14].real = 0.0;    v->c[14].imag = 0.0;
  v->c[15].real = 0.0;    v->c[15].imag = 0.0;
#if (DIMF > 16)
  v->c[16].real = 0.0;    v->c[16].imag = 0.0;
  v->c[17].real = 0.0;    v->c[17].imag = 0.0;
  v->c[18].real = 0.0;    v->c[18].imag = 0.0;
  v->c[19].real = 0.0;    v->c[19].imag = 0.0;
  v->c[20].real = 0.0;    v->c[20].imag = 0.0;
  v->c[21].real = 0.0;    v->c[21].imag = 0.0;
  v->c[22].real = 0.0;    v->c[22].imag = 0.0;
  v->c[23].real = 0.0;    v->c[23].imag = 0.0;
  v->c[24].real = 0.0;    v->c[24].imag = 0.0;
#if (DIMF > 25)
  v->c[25].real = 0.0;    v->c[25].imag = 0.0;
  v->c[26].real = 0.0;    v->c[26].imag = 0.0;
  v->c[27].real = 0.0;    v->c[27].imag = 0.0;
  v->c[28].real = 0.0;    v->c[28].imag = 0.0;
  v->c[29].real = 0.0;    v->c[29].imag = 0.0;
  v->c[30].real = 0.0;    v->c[30].imag = 0.0;
  v->c[31].real = 0.0;    v->c[31].imag = 0.0;
  v->c[32].real = 0.0;    v->c[32].imag = 0.0;
  v->c[33].real = 0.0;    v->c[33].imag = 0.0;
  v->c[34].real = 0.0;    v->c[34].imag = 0.0;
  v->c[35].real = 0.0;    v->c[35].imag = 0.0;
#if (DIMF > 36)
  v->c[36].real = 0.0;    v->c[36].imag = 0.0;
  v->c[37].real = 0.0;    v->c[37].imag = 0.0;
  v->c[38].real = 0.0;    v->c[38].imag = 0.0;
  v->c[39].real = 0.0;    v->c[39].imag = 0.0;
  v->c[40].real = 0.0;    v->c[40].imag = 0.0;
  v->c[41].real = 0.0;    v->c[41].imag = 0.0;
  v->c[42].real = 0.0;    v->c[42].imag = 0.0;
  v->c[43].real = 0.0;    v->c[43].imag = 0.0;
  v->c[44].real = 0.0;    v->c[44].imag = 0.0;
  v->c[45].real = 0.0;    v->c[45].imag = 0.0;
  v->c[46].real = 0.0;    v->c[46].imag = 0.0;
  v->c[47].real = 0.0;    v->c[47].imag = 0.0;
  v->c[48].real = 0.0;    v->c[48].imag = 0.0;
#if (DIMF > 49)
  register int i;
  for (i = 49; i < DIMF; i++) {
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;
  }
#endif
#endif
#endif
#endif
#endif
#endif
}
// -----------------------------------------------------------------
