/**
   @file
   @ingroup dsp
*/

#include "bpm_dsp.h"

int gaussian_filter_coeffs( filter_t *f ) {

  // implements the calculation of the gaussian filter coefficients
  // into the filter_t structure, the filter is implemented as a
  // non-causal FIR filter, in order to mimic the yury style filter
  // from the original ddc routines. Don't bother about the zero
  // calculation, just spit out the coefficients...
  int i = 0;
  double bw = 0.;
  char msg[80];

  if ( f->options & GAUSSIAN_SIGMA_BW ) {
    // f1 is frequency of the gaussian sigma bandwidth
    bw = f->f1; 
  } else {
    // f1 is frequency of -3dB bandwidth (default !!)
    bw = sqrt(SQR(f->f1)/(-2.*log(pow(10.,-3./20.)))); // convert bw to -3dB level
  }

  // number of causal coefficients
  f->nxc = (int) dround( sqrt( 2.*log( 2.*PI*bw / ( sqrt(2.*PI) * f->gauss_cutoff ) ) ) / 
			 (2.*PI*bw) * f->fs ) + 1 ;
  
  if ( ( f->nxc > MAXPZ ) || ( f->nxc >= f->ns ) ) {
    sprintf( msg, "Too many Gaussian coefficients : %d, encrease filter BW, or cutoff parameter",
	     f->nxc );
    bpm_error( msg, __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // number of anti-causal coefficients
  f->nxc_ac = f->nxc;
  

  // we kill some factors in the filter coefficient calculation since we
  // normalise afterwards anyway using the filter gain...
  for ( i = 0; i<f->nxc; i++ ) {
    // the causal coefficient array
    //f->xc[i] = 2.*PI*bw*exp(-SQR(2*PI*bw*(double) (f->nxc-1-i) /f->fs )/2.)/sqrt( 2.*PI );
    f->xc[i] = exp(-SQR(2*PI*bw*(double) (f->nxc-1-i) /f->fs )/2.);
    
    // the anti-causal coefficient array
    //f->xc_ac[i] = 2.*PI*bw*exp(-SQR(2*PI*bw*(double)( i )/f->fs)/2.)/sqrt( 2.*PI );
    f->xc_ac[i] = exp(-SQR(2*PI*bw*(double)( i )/f->fs)/2.);
  }
  
  // let's calculate the gain afterwards, we normalise anyway, so any factors in the 
  // filter coefficients don't really matter...
  f->gain = 0.;
  for ( i=0; i<f->nxc;i++)    f->gain += f->xc[i];
  for ( i=1; i<f->nxc_ac;i++) f->gain += f->xc_ac[i];

  return BPM_SUCCESS;
}
