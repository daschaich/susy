// -----------------------------------------------------------------
// Return the dot product of two vectors: adag b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/susy.h"

complex dot(vector *a, vector *b) {
#ifndef FAST
  register complex dot;
  dot.real = a->c[0].real * b->c[0].real + a->c[0].imag * b->c[0].imag;
  dot.imag = a->c[0].real * b->c[0].imag - a->c[0].imag * b->c[0].imag;

  dot.real += a->c[1].real * b->c[1].real + a->c[1].imag * b->c[1].imag;
  dot.imag += a->c[1].real * b->c[1].imag - a->c[1].imag * b->c[1].imag;

  dot.real += a->c[2].real * b->c[2].real + a->c[2].imag * b->c[2].imag;
  dot.imag += a->c[2].real * b->c[2].imag - a->c[2].imag * b->c[2].imag;

  dot.real += a->c[3].real * b->c[3].real + a->c[3].imag * b->c[3].imag;
  dot.imag += a->c[3].real * b->c[3].imag - a->c[3].imag * b->c[3].imag;
#if (DIMF > 4)
  dot.real += a->c[4].real * b->c[4].real + a->c[4].imag * b->c[4].imag;
  dot.imag += a->c[4].real * b->c[4].imag - a->c[4].imag * b->c[4].imag;

  dot.real += a->c[5].real * b->c[5].real + a->c[5].imag * b->c[5].imag;
  dot.imag += a->c[5].real * b->c[5].imag - a->c[5].imag * b->c[5].imag;

  dot.real += a->c[6].real * b->c[6].real + a->c[6].imag * b->c[6].imag;
  dot.imag += a->c[6].real * b->c[6].imag - a->c[6].imag * b->c[6].imag;

  dot.real += a->c[7].real * b->c[7].real + a->c[7].imag * b->c[7].imag;
  dot.imag += a->c[7].real * b->c[7].imag - a->c[7].imag * b->c[7].imag;

  dot.real += a->c[8].real * b->c[8].real + a->c[8].imag * b->c[8].imag;
  dot.imag += a->c[8].real * b->c[8].imag - a->c[8].imag * b->c[8].imag;
#if (DIMF > 9)
  dot.real += a->c[9].real * b->c[9].real + a->c[9].imag * b->c[9].imag;
  dot.imag += a->c[9].real * b->c[9].imag - a->c[9].imag * b->c[9].imag;

  dot.real += a->c[10].real * b->c[10].real + a->c[10].imag * b->c[10].imag;
  dot.imag += a->c[10].real * b->c[10].imag - a->c[10].imag * b->c[10].imag;

  dot.real += a->c[11].real * b->c[11].real + a->c[11].imag * b->c[11].imag;
  dot.imag += a->c[11].real * b->c[11].imag - a->c[11].imag * b->c[11].imag;

  dot.real += a->c[12].real * b->c[12].real + a->c[12].imag * b->c[12].imag;
  dot.imag += a->c[12].real * b->c[12].imag - a->c[12].imag * b->c[12].imag;

  dot.real += a->c[13].real * b->c[13].real + a->c[13].imag * b->c[13].imag;
  dot.imag += a->c[13].real * b->c[13].imag - a->c[13].imag * b->c[13].imag;

  dot.real += a->c[14].real * b->c[14].real + a->c[14].imag * b->c[14].imag;
  dot.imag += a->c[14].real * b->c[14].imag - a->c[14].imag * b->c[14].imag;

  dot.real += a->c[15].real * b->c[15].real + a->c[15].imag * b->c[15].imag;
  dot.imag += a->c[15].real * b->c[15].imag - a->c[15].imag * b->c[15].imag;
#if (DIMF > 16)
  dot.real += a->c[16].real * b->c[16].real + a->c[16].imag * b->c[16].imag;
  dot.imag += a->c[16].real * b->c[16].imag - a->c[16].imag * b->c[16].imag;

  dot.real += a->c[17].real * b->c[17].real + a->c[17].imag * b->c[17].imag;
  dot.imag += a->c[17].real * b->c[17].imag - a->c[17].imag * b->c[17].imag;

  dot.real += a->c[18].real * b->c[18].real + a->c[18].imag * b->c[18].imag;
  dot.imag += a->c[18].real * b->c[18].imag - a->c[18].imag * b->c[18].imag;

  dot.real += a->c[19].real * b->c[19].real + a->c[19].imag * b->c[19].imag;
  dot.imag += a->c[19].real * b->c[19].imag - a->c[19].imag * b->c[19].imag;

  dot.real += a->c[20].real * b->c[20].real + a->c[20].imag * b->c[20].imag;
  dot.imag += a->c[20].real * b->c[20].imag - a->c[20].imag * b->c[20].imag;

  dot.real += a->c[21].real * b->c[21].real + a->c[21].imag * b->c[21].imag;
  dot.imag += a->c[21].real * b->c[21].imag - a->c[21].imag * b->c[21].imag;

  dot.real += a->c[22].real * b->c[22].real + a->c[22].imag * b->c[22].imag;
  dot.imag += a->c[22].real * b->c[22].imag - a->c[22].imag * b->c[22].imag;

  dot.real += a->c[23].real * b->c[23].real + a->c[23].imag * b->c[23].imag;
  dot.imag += a->c[23].real * b->c[23].imag - a->c[23].imag * b->c[23].imag;

  dot.real += a->c[24].real * b->c[24].real + a->c[24].imag * b->c[24].imag;
  dot.imag += a->c[24].real * b->c[24].imag - a->c[24].imag * b->c[24].imag;
#if (DIMF > 25)
  dot.real += a->c[25].real * b->c[25].real + a->c[25].imag * b->c[25].imag;
  dot.imag += a->c[25].real * b->c[25].imag - a->c[25].imag * b->c[25].imag;

  dot.real += a->c[26].real * b->c[26].real + a->c[26].imag * b->c[26].imag;
  dot.imag += a->c[26].real * b->c[26].imag - a->c[26].imag * b->c[26].imag;

  dot.real += a->c[27].real * b->c[27].real + a->c[27].imag * b->c[27].imag;
  dot.imag += a->c[27].real * b->c[27].imag - a->c[27].imag * b->c[27].imag;

  dot.real += a->c[28].real * b->c[28].real + a->c[28].imag * b->c[28].imag;
  dot.imag += a->c[28].real * b->c[28].imag - a->c[28].imag * b->c[28].imag;

  dot.real += a->c[29].real * b->c[29].real + a->c[29].imag * b->c[29].imag;
  dot.imag += a->c[29].real * b->c[29].imag - a->c[29].imag * b->c[29].imag;

  dot.real += a->c[30].real * b->c[30].real + a->c[30].imag * b->c[30].imag;
  dot.imag += a->c[30].real * b->c[30].imag - a->c[30].imag * b->c[30].imag;

  dot.real += a->c[31].real * b->c[31].real + a->c[31].imag * b->c[31].imag;
  dot.imag += a->c[31].real * b->c[31].imag - a->c[31].imag * b->c[31].imag;

  dot.real += a->c[32].real * b->c[32].real + a->c[32].imag * b->c[32].imag;
  dot.imag += a->c[32].real * b->c[32].imag - a->c[32].imag * b->c[32].imag;

  dot.real += a->c[33].real * b->c[33].real + a->c[33].imag * b->c[33].imag;
  dot.imag += a->c[33].real * b->c[33].imag - a->c[33].imag * b->c[33].imag;

  dot.real += a->c[34].real * b->c[34].real + a->c[34].imag * b->c[34].imag;
  dot.imag += a->c[34].real * b->c[34].imag - a->c[34].imag * b->c[34].imag;

  dot.real += a->c[35].real * b->c[35].real + a->c[35].imag * b->c[35].imag;
  dot.imag += a->c[35].real * b->c[35].imag - a->c[35].imag * b->c[35].imag;
#if (DIMF > 36)
  dot.real += a->c[36].real * b->c[36].real + a->c[36].imag * b->c[36].imag;
  dot.imag += a->c[36].real * b->c[36].imag - a->c[36].imag * b->c[36].imag;

  dot.real += a->c[37].real * b->c[37].real + a->c[37].imag * b->c[37].imag;
  dot.imag += a->c[37].real * b->c[37].imag - a->c[37].imag * b->c[37].imag;

  dot.real += a->c[38].real * b->c[38].real + a->c[38].imag * b->c[38].imag;
  dot.imag += a->c[38].real * b->c[38].imag - a->c[38].imag * b->c[38].imag;

  dot.real += a->c[39].real * b->c[39].real + a->c[39].imag * b->c[39].imag;
  dot.imag += a->c[39].real * b->c[39].imag - a->c[39].imag * b->c[39].imag;

  dot.real += a->c[40].real * b->c[40].real + a->c[40].imag * b->c[40].imag;
  dot.imag += a->c[40].real * b->c[40].imag - a->c[40].imag * b->c[40].imag;

  dot.real += a->c[41].real * b->c[41].real + a->c[41].imag * b->c[41].imag;
  dot.imag += a->c[41].real * b->c[41].imag - a->c[41].imag * b->c[41].imag;

  dot.real += a->c[42].real * b->c[42].real + a->c[42].imag * b->c[42].imag;
  dot.imag += a->c[42].real * b->c[42].imag - a->c[42].imag * b->c[42].imag;

  dot.real += a->c[43].real * b->c[43].real + a->c[43].imag * b->c[43].imag;
  dot.imag += a->c[43].real * b->c[43].imag - a->c[43].imag * b->c[43].imag;

  dot.real += a->c[44].real * b->c[44].real + a->c[44].imag * b->c[44].imag;
  dot.imag += a->c[44].real * b->c[44].imag - a->c[44].imag * b->c[44].imag;

  dot.real += a->c[45].real * b->c[45].real + a->c[45].imag * b->c[45].imag;
  dot.imag += a->c[45].real * b->c[45].imag - a->c[45].imag * b->c[45].imag;

  dot.real += a->c[46].real * b->c[46].real + a->c[46].imag * b->c[46].imag;
  dot.imag += a->c[46].real * b->c[46].imag - a->c[46].imag * b->c[46].imag;

  dot.real += a->c[47].real * b->c[47].real + a->c[47].imag * b->c[47].imag;
  dot.imag += a->c[47].real * b->c[47].imag - a->c[47].imag * b->c[47].imag;

  dot.real += a->c[48].real * b->c[48].real + a->c[48].imag * b->c[48].imag;
  dot.imag += a->c[48].real * b->c[48].imag - a->c[48].imag * b->c[48].imag;
#if (DIMF > 49)
  register int i;
  for (i = 49; i < DIMF; i++) {
    dot.real += a->c[i].real * b->c[i].real + a->c[i].imag * b->c[i].imag;
    dot.imag += a->c[i].real * b->c[i].imag - a->c[i].imag * b->c[i].imag;
  }
#endif
#endif
#endif
#endif
#endif
#endif
  return dot;
#else // FAST version for DIMF=3 only
  register Real ar, ai, br, bi, cr, ci;
  register complex dot;

  ar = a->c[0].real;
  ai = a->c[0].imag;
  br = b->c[0].real;
  bi = b->c[0].imag;
  cr = ar * br + ai * bi;
  ci = ar * bi - ai * br;

  ar = a->c[1].real;
  ai = a->c[1].imag;
  br = b->c[1].real;
  bi = b->c[1].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  ar = a->c[2].real;
  ai = a->c[2].imag;
  br = b->c[2].real;
  bi = b->c[2].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  dot.real = cr;
  dot.imag = ci;
  return dot;
#endif
}
// -----------------------------------------------------------------
