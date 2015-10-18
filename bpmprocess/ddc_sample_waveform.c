/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <stdlib.h>

#include <bpm/bpm_units.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

#ifndef DDC_SAMPLE_DEBUG__
#define DDC_SAMPLE_DEBUG__
#endif

#undef DDC_SAMPLE_DEBUG__


int ddc_sample_waveform( doublewf_t *w, double frequency, filter_t *filt, int iSample,
			 double t0, double tdecay, double *amp, double *phase,
			 doublewf_t *buf_re, doublewf_t *buf_im ) {
  double tSample;

  if ( ! w ) {
    bpm_error( "Invalid waveform pointer in ddc_sample_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // calculate tSample from iSample

  bpm_error( "Not Yet Implemented... please use ddc_waveform instead !!!", __FILE__, __LINE__ );
  return BPM_FAILURE;


/* 
  // extrapolate back ....
  filt_ddc[Re] *= exp( ( tSample - t0 ) / tdecay );
  filt_ddc[Im] *= exp( ( tSample - t0 ) / tdecay );

  *amp   = sqrt( SQR( filt_ddc[Re] ) + SQR( filt_ddc[Im] ) );
  *phase = atan2( filt_ddc[Im], filt_ddc[Re] ); 
*/

  // Need routine apply_filter_at( filter, wf, t0+t0Offset );

  return 0;
}


/*
int ddc_sample_waveform( doublewf_t *w, int nbits, double t0, double t0Offset, 
			 double freq, double tdecay, double filtBW, double epsFilt,
			 double *amp, double *phase ) {
  
  int unsat_sample;
  double dt, tunsat, tstart, tfilter, tstop, omega;
  double filt_ddc[2];
  double **ddc, 
  char msg[80];
  int istart, ifilter, istop, i;

  if ( ! w ) {
    bpm_error( "Invalid waveform pointer in ddc_sample_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // get last unsaturated sample
  if ( check_saturation( w, nbits, &unsat_sample ) < 0 ) {
    bpm_error( "Error in check_saturation in ddc_sample_waveform(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

#ifdef DDC_SAMPLE_DEBUG__
  printf( "DBG:[ddc_sample_waveform] last unsaturated sample : %d\n", unsat_sample );
  printf( "DBG:[ddc_sample_waveform] t0Offset : %f\n", t0Offset );
  fflush( stdout );
#endif

  // allocate complex ddc waveform
  ddc = alloc_complex_wave_double( w->ns );
  if ( ! ddc ) {
    bpm_error( "Unable to allocate memory for ddc in ddc_sample_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // remove offset and downmix the waveform
  if ( downmix_waveform( w->wf, w->ns, w->fs, freq, ddc ) == BPM_FAILURE ) {
    bpm_error( "Error in downmixing in ddc_sample_waveform(...)", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // filter bandwidth given in MHz, convert it to rad/s
  omega = 2*PI*filtBW;
  
  // determine cut-off time from cut-off parameter epsFilt
  dt = sqrt( 2.*log( omega / ( sqrt(2.*PI) * epsFilt ) ) ) / omega;

  tstart  = t0Offset - dt;
  tfilter = t0Offset;
  tstop   = t0Offset + dt;

#ifdef DDC_SAMPLE_DEBUG__
  printf( "DBG:[ddc_sample_waveform] dt= %f, tstart= %f, tfilter= %f, tstop= %f\n", 
	  dt/usec, tstart/usec, tfilter/usec, tstop/usec );
  fflush( stdout );
#endif

  // set tstart according to last unsaturated sample
  //sample_to_time( fs, ns, unsat_sample, &tunsat ); 
  //if ( tunsat <= tstart ) { 
  //  time_to_sample( fs, ns, tstart, &istart );
  //} else {
  //  istart = unsat_sample;
  //}

  time_to_sample( fs, ns, tstart, &istart );
  time_to_sample( fs, ns, tfilter, &ifilter );
  time_to_sample( fs, ns, tstop, &istop );
  
  if (istart < unsat_sample) {
    istop += (unsat_sample - istart);
    ifilter += (unsat_sample - istart);
    istart = unsat_sample;
  }

#ifdef DDC_SAMPLE_DEBUG__
  printf( "DBG:[ddc_sample_waveform] istart= %d, ifilter= %d, istop= %d\n", 
	  istart, ifilter, istop );
  printf( "DBG:[ddc_sample_waveform] omega= %f, gamma= %f\n", freq, tdecay);
	  
  fflush( stdout );
#endif


  // we need only the value at tfilter, so just perform one gaussian filter step
  if ( ddc_gaussfilter_step( ddc, ns, fs, istart, istop, tfilter, filtBW, filt_ddc ) 
       == BPM_FAILURE ) {
    bpm_error( "Error in performing gaussian filter step in ddc_sample_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  filt_ddc[Re] *= exp( ( tfilter - t0 ) / tdecay );
  filt_ddc[Im] *= exp( ( tfilter - t0 ) / tdecay );

  *amp   = sqrt( SQR( filt_ddc[Re] ) + SQR( filt_ddc[Im] ) );
  *phase = atan2( filt_ddc[Im], filt_ddc[Re] ); 

  free_complex_wave_double( ddc, w->ns );

  return BPM_SUCCESS;
}
*/
