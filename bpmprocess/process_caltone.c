/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>
#include <bpm/bpm_dsp.h>

int process_caltone( doublewf_t *signal, bpmconf_t *bpm, bpmproc_t *proc, unsigned int mode ) {

  char msg[128];

  if ( ! bpm || ! signal || ! proc ) {
    bpm_error( "Invalid pointer arguments in process_caltone(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // do we have a signal ?
  if ( ! signal ) {
    sprintf( msg, "No signal present for BPM %s in process_waveform(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  /* ------------------------------- check for saturation ------------------------------------ */
  proc->saturated = check_saturation( signal, bpm->digi_nbits, &(proc->iunsat) );

  // report saturation in caltone... is not good :s
  if ( proc->saturated ) {
    bpm_warning( "Calibration tone is saturated, not updating caltone information...",
		 __FILE__, __LINE__ );
  } else {

    /* ------------------------------- subtract the pedestal ----------------------------------- */
    // determing voltage offset and amplitude noise in adc channels (pedestal)
    if ( get_pedestal( signal, 20, &(proc->voltageoffset) ,&(proc->ampnoise) ) == BPM_FAILURE ) {
      sprintf( msg, "Error getting pedestal of BPM %s in process_waveform(...)", bpm->name );
      bpm_error( msg, __FILE__, __LINE__ );
      return BPM_FAILURE;
    }
    
    // subtract the pedestal
    doublewf_bias( -proc->voltageoffset, signal );


    /* ------------------------------ check whether to do FFT ? -------------------------------- */
    if ( mode & PROC_DO_FFT ) {

      // compute the ft
      if ( fft_waveform( signal, proc->ft ) == BPM_FAILURE ) {
	sprintf( msg, "Could not perform fft for BPM %s in process_caltone(...)", bpm->name );
	bpm_warning( msg, __FILE__, __LINE__ );
	
      } else {
	proc->fft_success = TRUE;
	
	if ( mode & PROC_FIT_FFT ) {
	  // fft done, fit for frequency and tdecay
	  if( fit_fft( proc->ft, &(proc->fft_freq), &(proc->fft_tdecay), NULL, NULL ) 
	      == BPM_FAILURE ) {
	    sprintf( msg, "Could not fit the FFT for BPM %s in process_waveform(...)", bpm->name );
	    bpm_warning( msg, __FILE__, __LINE__ );
	  } 
	  
	}
	
      }

    } /* if ( mode & PROC_DO_FFT ) */


    /* ------------------------------ check whether to do DDC ? -------------------------------- */
    if ( mode & PROC_DO_DDC ) {
    
      // if we must to the full DDC, do it and sample afterwards
      if ( ddc_waveform( signal, bpm->ddc_ct_freq, bpm->ddc_ct_filter, proc->dc,
			 bpm->ddc_buffer_re, bpm->ddc_buffer_im ) == BPM_FAILURE ) {
	sprintf( msg, "Could not ddc BPM %s waveform in process_caltone(...)", bpm->name );
	bpm_warning( msg, __FILE__, __LINE__ );
      } else {
	proc->ddc_success = TRUE;
	
	// Now sample the ddc waveform at the iSample requested, extrapolate amplitude
 	proc->ddc_amp   = c_abs( proc->dc->wf[ bpm->ddc_ct_iSample ] );
	proc->ddc_phase = c_arg( proc->dc->wf[ bpm->ddc_ct_iSample ] );
	norm_phase( &(proc->ddc_phase) );

	// store phase and amplitude in the specially foreseen variables to store them
	// inbetween pulses...
	proc->ddc_ct_amp   = proc->ddc_amp;
	proc->ddc_ct_phase = proc->ddc_phase;
      }
      
    } /* if ( mode & PROC_DO_DDC ) */

  } // no saturation...

  return BPM_SUCCESS;
}
