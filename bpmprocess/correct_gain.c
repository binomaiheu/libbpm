/**
   @file
   @ingroup processing
*/

#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int correct_gain( bpmproc_t *proc, bpmcalib_t *cal, unsigned int mode ) {
  
  if ( ! proc || ! cal ) {
    bpm_error( "Invalid pointer arguments in correct_gain(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( proc->ddc_success ) {

    // scale the amplitude with the ratio of gain change between the caltone at the time
    // of calibration and the latest caltone amplitude in the bpmproc_t structure
    if ( ( mode & PROC_CORR_AMP ) && ( fabs( proc->ddc_ct_amp ) > 0. ) ) {
      proc->ddc_amp *= cal->ddc_ct_amp / proc->ddc_ct_amp;
    }
    // shift the phase by the phaseshift between the caltone at the time of calibration
    // and the lateste caltone phase in the bpmproc_t structure
    if ( mode & PROC_CORR_PHASE ) {
      proc->ddc_phase -= ( proc->ddc_ct_phase - cal->ddc_ct_phase );  
    }

  } /* if ( proc->ddc_success ) */


  if ( proc->fit_success ) {

    // same for the fitted stuf
    if ( ( mode & PROC_CORR_AMP ) && ( fabs( proc->fit_ct_amp ) > 0. ) ) {
      proc->fit_amp *= cal->fit_ct_amp / proc->fit_ct_amp;
    }

    if ( mode & PROC_CORR_PHASE ) {
      proc->fit_phase -= ( proc->fit_ct_phase - cal->fit_ct_phase );  
    }    
  } /* if ( proc->fit_success ) */

  return BPM_SUCCESS;
}
