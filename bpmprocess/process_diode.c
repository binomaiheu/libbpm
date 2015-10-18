/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int process_diode( doublewf_t *signal, bpmconf_t *conf, bpmproc_t *proc ) {

  static char msg[128];
  int i, ok = 0;
  wfstat_t s;
  
  if ( ! conf || ! signal || ! proc ) {
    bpm_error( "Invalid pointer arguments in process_diode(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! signal ) {
    sprintf( msg, "Invalid signal pointer for %s in process_diode(...)", conf->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( conf->cav_type == diode ) {
    // the signal supplied really is a diode pulse, so just fit it to determine
    // t0 !

    if ( fit_diodepulse( signal, &(proc->t0) ) == BPM_FAILURE ) {
      proc->t0 = 0.;
      bpm_error( "Could not fit the diode pulse in process_diode(...)", 
		 __FILE__, __LINE__ );
      return BPM_FAILURE;
    } 

  } else if ( conf->forced_trigger ) {
    // the signal supplied is a dipole or monopole pulse, that is abused
    // as trigger pulse !
    doublewf_basic_stats( signal, 0, 20, &s );
    for ( i=0; i<signal->ns; i++ ) {
      // we subtract the mean, since we don't know whether the signal was pedestal
      // subtracted...
      if ( fabs(signal->wf[i] - s.mean) > 10. * s.rms ) {
	proc->t0 = (double) i / signal->fs;
	ok = 1;
	break;
      }
    }
    // check it
    if ( ok != 1 ) {
      proc->t0 = 0.;
      sprintf( msg, "No onset of waveform found for %s in process_diode(...), pulse probably noise",
	       conf->name );
      bpm_error( msg, __FILE__, __LINE__ );
      return BPM_FAILURE;
    }
  } else {
    // this structure isn't a trigger signal or a forced trigger... are you really
    // sure you want to do this ???
    sprintf( msg, "Try to handle BPM %s through process_diode(...), don't think you want this...", 
	     conf->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  return BPM_SUCCESS;
}

