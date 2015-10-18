/**
   @file
   @ingroup dsp
*/
#include <bpm/bpm_dsp.h>


void norm_phase( double *phase ) {

  if ( ! phase ) return;
  while ( *phase <  0.0    ) *phase += 2.0*PI;
  while ( *phase >= 2.0*PI ) *phase -= 2.0*PI;

  return;
}
