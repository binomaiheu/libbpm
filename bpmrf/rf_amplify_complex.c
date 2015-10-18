/**
   @file
   @ingroup rf
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_nr.h>
#include <bpm/bpm_wf.h>

/**
   Amplifies the signal by the level dB. The voltage gain is calculated:
   @f[ gain = \sqrt{ 10^{ \frac{db}{20} } } @f]
   @param RF waveform to be processed
   @param dB gain (or attenuation) in dB
   @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
*/
int rf_amplify_complex( complexwf_t *RF, double dB ) {

  complex_t f;
 
  if ( ! RF ) {
    bpm_error( "Invalid pointer arguments in rf_amplify(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  f.re = pow ( 10., dB / 20. );
  f.im = 0.;

  complexwf_scale( f, RF );

  return BPM_SUCCESS;
}
// end of file
