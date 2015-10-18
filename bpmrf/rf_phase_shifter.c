/**
   @file
   @ingroup rf
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_nr.h>
#include <bpm/bpm_wf.h>

/**
   Rotates the phase of the signal by the amount specified
   @param RF waveform to be processed
   @param rotation phase rotation in degrees
   @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
*/
int rf_phase_shifter( complexwf_t *RF, double rotation ) {

  double phi;
  complex_t f;
 
  if ( ! RF ) {
    bpm_error( "Invalid pointer arguments in rf_amplify(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  phi = PI * rotation/180.;
  f.re = cos (phi);
  f.im = sin (phi);

  complexwf_scale( f, RF );

  return BPM_SUCCESS;
}
// end of file
