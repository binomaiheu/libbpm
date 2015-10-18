/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>


int process_monopole( doublewf_t *sig, bpmconf_t *bpm, bpmcalib_t *cal, 
		      bpmproc_t *proc, bpmproc_t *trig, unsigned int mode ) {

  char msg[128];
 
  if ( ! sig || ! bpm || ! cal || ! proc || ! trig ) {
    bpm_error( "Invalid pointer arguments in process_monopole(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( process_waveform( sig, bpm, proc, trig, mode ) == BPM_FAILURE ) {
    sprintf( msg, "Unable to process waveform for BPM %s in process_monopole(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;

  }

  if ( correct_gain( proc, cal, mode ) == BPM_FAILURE ) {
    sprintf( msg, "Unable to correct gains for BPM %s in process_monopole(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// end of file
