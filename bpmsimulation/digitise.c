/**
   @file
   @ingroup sim
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_simulation.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_nr.h>
#include <bpm/bpm_wf.h>

int digitise( doublewf_t *IF, int nbits, double range_min, double range_max,
	      double clock_jitter, double digi_noise, unsigned int ipmode, 
	      intwf_t *wf ) {
  
  double range, ti, amp;
  int    i;

  if (( ! IF ) || ( ! wf )) {
    bpm_error( "Invalid pointer arguments in digitise(...)", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  if ( nbits <= 0 ) {
    bpm_error( "Invalid number of ADC bits in digitise(...)",
	       __FILE__, __LINE__);
    return BPM_FAILURE;
  }

  if ( range_max <= range_min ) {
    bpm_error( "Invalid range setting in digitise(...)",
               __FILE__, __LINE__);
    return BPM_FAILURE;
  }

  range = pow( 2., nbits ); // range of the ADC
  
  // now at first sample
  for ( i = 0; i < wf->ns; i++ ) {

    // ADC timestamp of this sample, (without time jitter)
    ti = (double) i / wf->fs;
    
    // add time jitter here when required
    if ( clock_jitter != 0 ) ti = nr_rangauss( ti, clock_jitter );

    // are we already beyond the IF time ?
    if ( ti > (double) (IF->ns-1) / IF->fs ) {
      // We are already past the last sample in the IF waveform
      // In this case we just fill the array with the pedestal 
      wf->wf[i] = (int) range / 2;
    } else {

      // resample the waveform
      amp = doublewf_getvalue( IF, ti, ipmode );

      // now convert this value to an ADC value, adding half of the range
      // as pedestal
      wf->wf[i] = (int) ( amp * range / ( range_max - range_min ) + range / 2 );
    }

  }

  // Add digitiser noise
  intwf_add_ampnoise( wf, digi_noise );

  // Address saturation in the ADC 
  for ( i = 0; i < wf->ns; i++ ) {
    if( wf->wf[i] < 0 ) wf->wf[i] = 0;
    if( wf->wf[i] > ( range - 1 ) ) wf->wf[i] = range - 1;
  }

  return BPM_SUCCESS;
}

// end of file
