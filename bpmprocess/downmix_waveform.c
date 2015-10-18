/**
   @file
   @ingroup processing
*/

#include <math.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int downmix_waveform( doublewf_t *w, double freq, complexwf_t *out ) {

  int i;

  if ( ! w || ! out ) {
    bpm_error( "Invalid pointer arguments in downmix_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // dowmixing
  for ( i=0; i<w->ns; i++ ) {
    out->wf[i].re = w->wf[i] * cos( 2. * PI * freq * (double) i / w->fs );
    out->wf[i].im = w->wf[i] * sin( 2. * PI * freq * (double) i / w->fs );
  }

  return BPM_SUCCESS;
}
