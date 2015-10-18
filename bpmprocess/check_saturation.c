/**
   @file
   @ingroup processing
*/
#include <math.h>
#include <limits.h>

#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int check_saturation( doublewf_t *w, int nbits, int *iunsat ) {

  int i;
  double min_value, max_value;
  double thresh = 15.;

  *iunsat = -INT_MAX;

  if ( ! w ) {
    bpm_error( "Invalid waveform pointer in check_saturation(...)",
	       __FILE__, __LINE__ );
    return -1;
  }

  min_value = thresh;
  max_value = (double) ( 1 << nbits ) - thresh;

  if ( max_value <= min_value ) {
    bpm_error( "Check number of bits in ADC and threshold for check_saturation(...)",
	       __FILE__, __LINE__ );
    return -1;
  }

  // search for sample
  for ( i = w->ns-1; i >= 0; --i ) {
    if ( ( w->wf[i] > max_value ) || ( w->wf[i] < min_value ) ) break;
  }

  // found saturation
  if ( i > 0 ) {
    if ( i < w->ns - 1 ) *iunsat = i + 1;
    return 1;
  }

  // no saturation
  *iunsat = 0;
  return 0;
}
