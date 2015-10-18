/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>
#include <bpm/bpm_dsp.h>


int process_waveform( doublewf_t *signal, bpmconf_t *bpm, bpmproc_t *proc, bpmproc_t *trig, 
		      unsigned int mode ) {
  
  int have_fft_pars = FALSE;
  int have_fit_pars = FALSE;
  double frequency, tdecay, amp, phase;
  char msg[128];
  double t0, tSample;
  
  if ( ! bpm || ! signal || ! proc ) {
    bpm_error( "Invalid pointer arguments in process_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // reset the success values for this pulse...
  proc->fft_success = FALSE;
  proc->fit_success = FALSE;
  proc->ddc_success = FALSE;

  // do we have a signal ?
  if ( ! signal ) {
    sprintf( msg, "No signal present for BPM %s in process_waveform(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  /* ------------------------------- check for saturation ------------------------------------ */
  proc->saturated = check_saturation( signal, bpm->digi_nbits, &(proc->iunsat) );


  /* ------------------------------- subtract the pedestal ----------------------------------- */
  // determing voltage offset and amplitude noise in adc channels (pedestal)
  if ( get_pedestal( signal, 20, &(proc->voltageoffset) ,&(proc->ampnoise) ) == BPM_FAILURE ) {
    sprintf( msg, "Error getting pedestal of BPM %s in process_waveform(...)", bpm->name );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
    
  // subtract the pedestal
  doublewf_bias( -proc->voltageoffset, signal );


  /* ----------------------------------- set the t0 ------------------------------------------ */
  if ( trig ) {
    // internal clock, we have a trigger
    proc->t0 = trig->t0;
  } else {
    // external clock, no trigger specified, assume t0 is set in configuration
    proc->t0 = bpm->t0;
  }

  /* ------------------------------ check whether to do FFT ? -------------------------------- */
  if ( mode & PROC_DO_FFT ) {

    // compute the ft
    if ( fft_waveform( signal, proc->ft ) == BPM_FAILURE ) {
      sprintf( msg, "Could not perform fft for BPM %s in process_waveform(...)", bpm->name );
      bpm_warning( msg, __FILE__, __LINE__ );

    } else {
      proc->fft_success = TRUE;
      
      if ( mode & PROC_FIT_FFT ) {
	// fft done, fit for frequency and tdecay
	if( fit_fft( proc->ft, &(proc->fft_freq), &(proc->fft_tdecay), NULL, NULL ) 
	    == BPM_FAILURE ) {
	  sprintf( msg, "Could not fit the FFT for BPM %s in process_waveform(...)", bpm->name );
	  bpm_warning( msg, __FILE__, __LINE__ );
	} else { 
	  have_fft_pars     = TRUE;
	}

      }

    }

  } /* if ( mode & PROC_DO_FFT ) */


  /* ------------------------------ check whether to do FIT ? -------------------------------- */
  if ( mode & PROC_DO_FIT ) {
    
    /*
      This stuff is a bit shaky... need to clean & test in more detail if needed,
      rely on ROOT/Minuit based fitting for this !!!
    */
    bpm_warning( "Libbpm waveform fitting is highly EXPERIMENTAL and contains bugs!",
		 __FILE__, __LINE__ );

    // initial parameters
    frequency = tdecay = amp = phase = 0.;

    // Initial parameter for the frequency ::
    if ( proc->fit_freq <= 0.* _MHz__ || proc->fit_freq > ( bpm->digi_freq/2.) ) {
      // BM. 11.02.2008
      // cut this out, don't want it to depend on the rf_LOfreq just yet...
      /*
	if ( ( bpm->rf_LOfreq > 0. ) && ( bpm->cav_freq > 0. ) ) {
	freq = ABS( bpm->cav_freq - bpm->rf_LOfreq );
	} else 
      */
      frequency = 20.0 * _MHz__;
    } else frequency = proc->fit_freq; // use last fitted value

    // Initial parameter for the decay constant
    if ( proc->fit_tdecay <= 10.* _nsec__ || proc->fit_tdecay > 1.0 * _sec__ ) {
      if ( bpm->cav_decaytime > 0. ) tdecay = bpm->cav_decaytime;
      else tdecay = 0.2 * _sec__;
    } else tdecay = proc->fit_tdecay; // use last fitted value

    // probably can think of something more clever here...
    if( ABS( proc->fit_amp ) < 100000. ) amp = proc->fit_amp; else amp = 3000.;
    phase  = proc->fit_phase;
    
    // if we have down the FFT and got the frequency from there, use that as an inital parameter
    // so overwrite the previous
    if( have_fft_pars ) {
      frequency   = proc->fft_freq;
      tdecay      = proc->fft_tdecay; 
    }

    // Slight hack to sort out saturation and trigger offsets
    // Change t0, fit, then extrapolate back in both amplitude and phase
    sample_to_time( bpm->digi_freq, bpm->digi_nsamples, proc->iunsat, &t0 );
    t0 = ( t0 > ( proc->t0 + bpm->fit_tOffset ) ? t0 : proc->t0 + bpm->fit_tOffset);

    if ( fit_waveform( signal, t0, frequency, tdecay, amp, phase,
		       &(proc->fit_freq), &(proc->fit_tdecay), &(proc->fit_amp) ,
		       &(proc->fit_phase) ) == BPM_FAILURE ) {
	sprintf( msg, "Could not fit for BPM %s in process_waveform(...)", bpm->name );
	bpm_warning( msg, __FILE__, __LINE__ );
    } else { 
      have_fit_pars = TRUE;
      proc->fit_success = TRUE;

      // Extrapolate back to the real t0 (proc->t0) !!
      proc->fit_amp   *= exp( ( t0 - proc->t0 ) / proc->fit_tdecay );
      proc->fit_phase -= (t0 - proc->t0) * proc->fit_freq * 2. * PI;
      norm_phase( &(proc->fit_phase) );
    }

  } /* if ( mode & PROC_DO_FIT ) */


  /* ------------------------------ check whether to do DDC ? -------------------------------- */
  if ( mode & PROC_DO_DDC ) {
    
    // initial parameters
    frequency = tdecay = 0.;

    // where does the ddc get its freq and tdecay parameters from ?
    if ( mode & PROC_DDC_FITFREQ ) {
      if ( have_fit_pars ) frequency = proc->fit_freq;
      else frequency = bpm->ddc_freq; // asked for fit but it failed :(
    } else if ( mode & PROC_DDC_FFTFREQ ) {
      if ( have_fft_pars ) frequency = proc->fft_freq;
      else frequency = bpm->ddc_freq; // asked for fft but it failed :(
    } else {
      // use default
      frequency = bpm->ddc_freq;
    }

    // same for tdecay
    if ( mode & PROC_DDC_FITTDECAY ) {
      if ( have_fit_pars ) tdecay = proc->fit_tdecay;
      else tdecay = bpm->ddc_tdecay; // asked for fit but it failed :(
    } else if ( mode & PROC_DDC_FFTTDECAY ) {
      if ( have_fft_pars ) tdecay = proc->fft_tdecay;
      else tdecay = bpm->ddc_tdecay; // asked for fft but it failed :(
    } else {
      // use default
      tdecay = bpm->ddc_tdecay;
    }
    
    // determine the sample time, handle saturation here...
    if ( proc->saturated == 1 ) {
      // We have saturation, get the index of the last unsaturated sample
      sample_to_time( bpm->digi_freq, bpm->digi_nsamples, proc->iunsat, &(proc->ddc_tSample) );
      // normally we should add here something which takes into account the bandwidth
      // of the filter since we shift the filter window, let's just take 
    } else {
      // The normall case, the t0 (trigger or fixed + the offset)
      proc->ddc_tSample = proc->t0 + bpm->ddc_tOffset; 
    }

    // convert to sample index !
    time_to_sample( bpm->digi_freq, bpm->digi_nsamples, proc->ddc_tSample, &(proc->ddc_iSample) );

    if ( mode & PROC_DDC_FULL ) {
      // if we must to the full DDC, do it and sample afterwards
      if ( ddc_waveform( signal, frequency, bpm->ddc_filter, proc->dc,
			 bpm->ddc_buffer_re, bpm->ddc_buffer_im ) == BPM_FAILURE ) {
	sprintf( msg, "Could not ddc BPM %s waveform in process_waveform(...)", bpm->name );
	bpm_warning( msg, __FILE__, __LINE__ );
	proc->ddc_amp   = 0.;
	proc->ddc_phase = 0.;

      } else {
	proc->ddc_success = TRUE;

	// Now sample the ddc waveform at the iSample requested, extrapolate amplitude
	proc->ddc_amp   = exp( ( proc->ddc_tSample - proc->t0 ) / tdecay ) * 
	  c_abs( proc->dc->wf[ proc->ddc_iSample ] );
	proc->ddc_phase = c_arg( proc->dc->wf[ proc->ddc_iSample ] );
	norm_phase( &(proc->ddc_phase) );
      }

    } else {
      // Execute the "short" ddc algorithm to save cpu cycles
      if ( ddc_sample_waveform( signal, frequency, bpm->ddc_filter, proc->ddc_iSample, 
				proc->t0, tdecay, &(proc->ddc_amp), &(proc->ddc_phase),
				bpm->ddc_buffer_re, bpm->ddc_buffer_im )  == BPM_FAILURE ) {
	sprintf( msg, "Could not ddc_sample BPM %s waveform in process_waveform(...)", bpm->name );
	bpm_warning( msg, __FILE__, __LINE__ );
      } else {
	proc->ddc_success = TRUE;
	norm_phase( &(proc->ddc_phase) );
      }
    }

  } /* if ( mode & PROC_DO_DDC ) */

  return BPM_SUCCESS;
}
