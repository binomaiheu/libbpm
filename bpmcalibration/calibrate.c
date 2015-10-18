/**
   @file
   @ingroup calib
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_calibration.h>
#include <bpm/bpm_nr.h>

int calibrate( bpmconf_t *bpm, bunchconf_t *bunch, bpmproc_t *proc, int npulses, 
	       bpmcalib_t *cal ) {

  if ( ! cal || ! proc || ! bunch ) {
    bpm_error( "Invalid pointer arguments in calibrate(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // BM: 07-04-2008
  // this needs to be re-written....

/*

  // First, fit a staright line to the data to get the IQ phase
  double *ival = calloc( npulses, sizeof(double) );
  double *qval = calloc( npulses, sizeof(double) );
  int i;

  for (i = 0; i < npulses; i++) {
    ival[i] = proc[i].ddc_I;
    qval[i] = proc[i].ddc_Q;
  }

  double intercept, gradient, siga, sigb, chi2, q;
  if ( nr_fit( ival, qval, npulses, NULL, 0, &intercept, &gradient, 
	       &siga, &sigb, &chi2, &q ) == BPM_FAILURE ) {
    bpm_error( "Straight line fit failed in calibrate(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  cal->IQphase = gradient;

  // Now calculate the scale
  double pos;
  double *scale = calloc( npulses, sizeof(double) );

  for (i = 0; i < npulses; i++) {

    get_pos( proc[i].ddc_Q, proc[i].ddc_I, cal->IQphase, 1., &pos );

    if ( bpm->cav_polarisation == horiz ) {
      scale[i] = bunch[i].bpmposition[0] / pos;
    } else {
      scale[i] = bunch[i].bpmposition[1] / pos;
    }
 
  }

  cal->posscale = nr_median(npulses, scale);
*/
  return BPM_SUCCESS;
}


