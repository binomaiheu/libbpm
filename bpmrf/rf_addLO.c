/**
   @file
   @ingroup rf
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_nr.h>
#include <math.h>
#include <bpm/bpm_wf.h>

/**
  Generates an LO waveform
  @param amp amplitude of the LO signal in Volts
  @param lofreq LO frequency
  @type locked or freerunning oscillator
  @phase phase of the signal, ignored if type is not "locked"
  @phasenoise phase noise to be added to the waveform
  @param LO generated waveform
  @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
*/
int rf_addLO( double amp, double lofreq, enum bpmphase_t type,
              double phase, double phasenoise, doublewf_t *LO ) {

  if ( ! LO ) {
    bpm_error( "Invalid LO pointer in rf_addLO(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // Generate the LO upon the phase type
  
  if ( type == locked ) {

    doublewf_add_cwtone( LO, amp, phase, lofreq, phasenoise);

  } else {

    doublewf_add_cwtone( LO, amp, nr_ranuniform( 0., 2*PI ), lofreq, phasenoise);

  }

  return BPM_SUCCESS;
}
// end of file
