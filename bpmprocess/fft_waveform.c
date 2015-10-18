/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <stdlib.h>

#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>
#include <bpm/bpm_dsp.h>

int fft_waveform( doublewf_t *w, complexwf_t *fft ) {

  if ( ! w || ! fft ) {
    bpm_error( "Invalid waveform pointers in fft_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return realfft( w, FFT_FORWARD, fft );
}

