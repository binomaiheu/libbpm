/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <stdlib.h>

#include <bpm/bpm_units.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int ddc_waveform( doublewf_t *w, double frequency, filter_t *filt, complexwf_t *dc,
		  doublewf_t *buf_re, doublewf_t *buf_im  ) {

  if ( ! w || ! filt || ! dc ) {
    bpm_error( "Invalid waveform pointer in ddc_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  return ddc( w, frequency, filt, dc, buf_re, buf_im );
}
