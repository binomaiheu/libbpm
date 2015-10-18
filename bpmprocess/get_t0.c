/**
   @file
   @ingroup processing
   Declared two helper routines which find the start and end samples for the fit...
*/
#include <stdlib.h>
#include <math.h>

#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>
#include <bpm/bpm_nr.h>

void _find_t0_startfit( double *wf, double ped, 
			int peak_sample, double peak_value, double peak_fraction, 
			int *start_sample ) {
  // starts from the peak sample and goes back in the waveform as long as the
  // value stays > the given fraction of the peak value
  int S = 0;
  for ( S = peak_sample; S > 0; S-- ) {
    if ( ( ABS( wf[S] - ped ) > peak_fraction * peak_value ) && 
	 ( ABS( wf[S] - ped ) > 6.5 ) ) {
      *start_sample = S;
    }
  }
  return;
} /* find_t0_startfit(...) */

// ---------------------------------------------------------------------------
  
void _find_t0_endfit( double *wf, double ped,
		      int peak_sample, double peak_value, double peak_fraction,
		      int *end_sample ) {
  // starts from 0 to the peak sample and advances in the waveform untill
  // the value is the given fraction of the peak value
  int S = 0;
  for ( S = 0; S < peak_sample; S++ ) {
    if ( ( ABS( wf[S] - ped ) ) < peak_fraction*peak_value ) {
      *end_sample = S;
    }
  }
  return;
} /* find_t0_endfit(...) */

// ---------------------------------------------------------------------------

int get_t0( doublewf_t *signal, double *t0 ) {
  

  double adc_ped, adc_err, peak_value, sample_value, start_value;
  double gradient, intercept, *xval, *yval;
  int peak_sample, start_sample, counter, end_sample, i, retcode;
  double siga, sigb, chi2, q;

  *t0 = -DBL_MAX; // some initial value :)

  if ( ! signal || ! t0 ) {
    bpm_error( "Invalid pointer arguments in get_t0(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // Find some useful info
  if( get_pedestal( signal, 20, &adc_ped, &adc_err ) == BPM_FAILURE ) {
    bpm_error( "Unable to retreive pedestal in get_t0(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // Find the peak sample
  // Note - cut off the last 10 samples in case there is a problem 
  //        with the end of the waveform
  peak_value = adc_err * 4.;
  peak_sample = 0;

  for ( i = 0; i < signal->ns - 10; ++i ) {
    sample_value = ABS( signal->wf[i] - adc_ped );
    
    if ( sample_value > peak_value ) {
      peak_value  = sample_value;
      peak_sample = i;
    }

  } /* for ( sample = 0; sample < ns - 10; ++sample ) */

  // Did we find a peak?
  if ( peak_sample == 0 ) { // No!
    bpm_error( "Could not find a peak in get_t0(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // Find the start point
  _find_t0_endfit( signal->wf, adc_ped, peak_sample, peak_value, 0.9, &end_sample );
  _find_t0_startfit( signal->wf, adc_ped, peak_sample, peak_value, 0.1, &start_sample );
  
  // Only one sample so try to find more
  if ( start_sample == end_sample ) {
    if ( bpm_verbose ) bpm_warning( "First fit initialisation failed, trying second",
				    __FILE__, __LINE__ );
    _find_t0_endfit( signal->wf, adc_ped, peak_sample, peak_value, 0.95, &end_sample );
    _find_t0_startfit( signal->wf, adc_ped, peak_sample, peak_value, 0.05, &start_sample );
  } 

  // Still only one point? Find more!
  if ( start_sample == end_sample ) {
    if ( bpm_verbose ) bpm_warning( "Second fit initialisation failed, trying third",
				    __FILE__, __LINE__ );
    _find_t0_endfit( signal->wf, adc_ped, peak_sample, peak_value, 0.975, &end_sample );
    _find_t0_startfit( signal->wf, adc_ped, peak_sample, peak_value, 0.025, &start_sample );
  } 

  // Still only one point? Give up...
  if ( ( end_sample == start_sample ) || 
       ( start_sample > end_sample  ) ) { // No!
    bpm_warning( "Cannot initialise fit, returning end_sample in get_t0(...)",
                 __FILE__, __LINE__ );
    *t0 = (double) end_sample / signal->fs; // return in time units
    return BPM_SUCCESS;
  }
 
  // Fit the intervening samples to a straight line and
  // return the fraction point
  xval = (double*) calloc( end_sample - start_sample + 1, sizeof(double) );
  yval = (double*) calloc( end_sample - start_sample + 1, sizeof(double) );

  if ( ! ( xval && yval ) ) {
    bpm_error( "Coudn't allocate memory in get_t0(...)", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  counter = 0;
  for ( i = start_sample; i < end_sample + 1; i++ )  {
    yval[counter]   = signal->wf[i] - adc_ped;
    xval[counter++] = (double) i;
  }

  // perform the fit
  if ( nr_fit( xval, yval, end_sample - start_sample + 1, NULL, 0, &intercept, &gradient, 
	       &siga, &sigb, &chi2, &q ) == BPM_FAILURE ) {
    bpm_error( "T0 straight line fit failed in get_t0(...)", __FILE__, __LINE__ );
    *t0 = -DBL_MAX;
    retcode = BPM_FAILURE;
  } else {
    // fit was successfull
    if( ABS( gradient ) > 0. ) {
      *t0 = ( peak_value * 0.5 - intercept ) / gradient;
      if ( *t0 < 0 || *t0 > signal->ns ) {
	bpm_error( "Fitted t0 value out of range in get_t0(...)",
		   __FILE__, __LINE__ );
	*t0 = -DBL_MAX;
	retcode = BPM_FAILURE;
      } else {
	// we have a good and proper t0, convert to time units...
	*t0 /= signal->fs; 
	retcode = BPM_SUCCESS;
      }
    } else {
      bpm_error( "Gradient in t0 fit equals 0 in get_t0(...)",
                 __FILE__, __LINE__ );
      *t0 = -DBL_MAX;
      retcode = BPM_FAILURE;
    }

  }
  
  // free the memory
  free( xval );
  free( yval );

  return retcode;
}
