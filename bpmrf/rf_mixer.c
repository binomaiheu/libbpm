/**
   @file
   @ingroup rf
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_wf.h>

/**
  Simulates an ideal mixer
  @param RF signal to mix
  @param LO local oscillator signal to mix with
  @param IF resulting signal containing the up and down
  converted terms
  @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
*/
int rf_mixer( doublewf_t *RF, doublewf_t *LO, doublewf_t *IF ) {


  if ( ! RF || ! LO || ! IF ) {
    bpm_error( "Invalid pointer arguments in rf_mixer(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // make copy of RF into the IF waveform, RF gets preserved
  doublewf_copy( IF, RF );

  // mix down
  doublewf_multiply( IF, LO );

  return BPM_SUCCESS;
}
// end of file
