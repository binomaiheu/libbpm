/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>


int process_dipole( doublewf_t *sig, bpmconf_t *bpm, bpmcalib_t *cal, bpmproc_t *proc,
		    bpmproc_t *trig, bpmproc_t *ampref, bpmproc_t *phaseref, unsigned int mode ) {
  
  char msg[128];

  if ( ! bpm || ! sig || ! cal || ! proc || ! trig || ! ampref || ! phaseref ) {
    bpm_error( "Invalid pointer arguments in process_dipole(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
    
  if ( process_waveform( sig, bpm, proc, trig, mode ) == BPM_FAILURE ) {
    sprintf( msg, "Unable to process waveform for BPM %s in process_dipole(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // correct gain according to the mode...
  if ( correct_gain( proc, cal, mode ) == BPM_FAILURE ) {
    sprintf( msg, "Unable to correct gains for BPM %s in process_dipole(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( postprocess_waveform( bpm, proc, cal, ampref, phaseref, mode ) == BPM_FAILURE ) {
    sprintf( msg, "Unable to handle post processing for BPM %s in process_dipole(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// end of file
