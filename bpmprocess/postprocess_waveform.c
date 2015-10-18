/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>
#include <bpm/bpm_dsp.h>

int postprocess_waveform(  bpmconf_t *bpm, bpmproc_t *proc, bpmcalib_t *cal, 
			   bpmproc_t *ampref, bpmproc_t *phaseref, unsigned int mode ) {
  char msg[80];
  
  if ( ! bpm || ! proc || ! cal || ! ampref || ! phaseref ) {
    bpm_error( "Invalid pointer arguments in postprocess_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // Initialise
  proc->ddc_I = 0.;
  proc->ddc_Q = 0.;
  proc->ddc_pos = 0.;
  proc->ddc_slope = 0.;

  if ( proc->ddc_success && ampref->ddc_success && phaseref->ddc_success ) {

    // Process info for DDC
    if ( get_IQ( proc->ddc_amp, proc->ddc_phase, ampref->ddc_amp, phaseref->ddc_phase, 
                 &(proc->ddc_Q), &(proc->ddc_I) ) == BPM_FAILURE ) {
      sprintf( msg, "Could not calculate I and Q for BPM %s in process_dipole(...)", bpm->name );
      bpm_error( msg, __FILE__, __LINE__ );
      proc->ddc_Q = 0.;
      proc->ddc_I = 0.;
      return BPM_FAILURE;
    } 

    // Store the corrected phase if required
    if ( !( mode & PROC_RAW_PHASE ) ) {
      proc->ddc_phase -= phaseref->ddc_phase;
      norm_phase( &(proc->ddc_phase) );
    }

    // Calculate position and slope from I & Q
    if ( get_pos( proc->ddc_Q, proc->ddc_I, cal->ddc_IQphase, cal->ddc_posscale, 
                  &(proc->ddc_pos) ) == BPM_FAILURE ) {
      sprintf( msg, "Could not get ddc position for BPM %s in process_dipole(...)", bpm->name );
      bpm_error( msg, __FILE__, __LINE__ );
      proc->ddc_pos = 0.;
      return BPM_FAILURE;
    }
    if ( get_slope( proc->ddc_Q, proc->ddc_I, cal->ddc_IQphase, cal->ddc_slopescale, 
                    &(proc->ddc_slope) ) == BPM_FAILURE ) {
      sprintf( msg, "Could not get ddc slope for BPM %s in process_dipole(...)", bpm->name );
      proc->ddc_slope = 0.;
      bpm_error( msg, __FILE__, __LINE__ );
      return BPM_FAILURE;
    }

  } /* if ( proc->ddc_success ) */


  // Same but for the fitted amplitude and phase...
  proc->fit_I = 0.;
  proc->fit_Q = 0.;
  proc->fit_pos = 0.;
  proc->fit_slope = 0.;

  if ( proc->fit_success && ampref->fit_success && phaseref->fit_success ) {

    // Process info for FIT
    if ( get_IQ( proc->fit_amp, proc->fit_phase, ampref->fit_amp, phaseref->fit_phase, 
                 &(proc->fit_Q), &(proc->fit_I) ) == BPM_FAILURE ) {
      sprintf( msg, "Could not calculate I and Q for BPM %s in process_dipole(...)", bpm->name );
      bpm_error( msg, __FILE__, __LINE__ );
      proc->fit_Q = 0.;
      proc->fit_I = 0.;
      return BPM_FAILURE;
    } 

    // Store the corrected phase if required
    if ( !( mode & PROC_RAW_PHASE ) ) {
      proc->fit_phase -= phaseref->fit_phase;
      norm_phase( &(proc->fit_phase) );
    }

    // Calculate position and slope from I & Q
    if ( get_pos( proc->fit_Q, proc->fit_I, cal->fit_IQphase, cal->fit_posscale, 
                  &(proc->fit_pos) ) == BPM_FAILURE ) {
      sprintf( msg, "Could not get fit position for BPM %s in process_dipole(...)", bpm->name );
      bpm_error( msg, __FILE__, __LINE__ );
      proc->fit_pos = 0.;
      return BPM_FAILURE;
    }
    if ( get_slope( proc->fit_Q, proc->fit_I, cal->fit_IQphase, cal->fit_slopescale, 
                    &(proc->fit_slope) ) == BPM_FAILURE ) {
      sprintf( msg, "Could not get fit slope for BPM %s in process_dipole(...)", bpm->name );
      proc->fit_slope = 0.;
      bpm_error( msg, __FILE__, __LINE__ );
      return BPM_FAILURE;
    }

  } /* if ( proc->fit_success ) */

  return BPM_SUCCESS;
}
  
