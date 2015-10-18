/**
   @file
   @ingroup rf
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_units.h>

/**
   Rectifies the given waveform assuming a single diode
   @param D the rectified signal
   @param RF the complex waveform to rectify
   @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
*/
int rf_rectify( doublewf_t *D, complexwf_t *RF ) {

  
  if ( ( ! RF) || ( ! D) ) {
    bpm_error( "Invalid IF pointer in rf_rectify(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  complexwf_getreal( D, RF);

  int i;

  for ( i = 0; i < D->ns; i++ ) {

    if ( D->wf[i] < 0 ) D->wf[i] = 0.;

  }  return BPM_SUCCESS;
}
// end of file
