#include "bpm/bpm_dsp.h"


void _shift_down( double *x, int n ) {
  int k=0;
  for( k=0; k<n-1; k++ ) x[k] = x[k+1];
  return;
}

// -------------------------------------------------------------------

void _reset( double *x ) {
  int i = 0;
  for( i=0;i<=MAXPZ;i++) x[i] = 0.;
}

// -------------------------------------------------------------------

int apply_filter( filter_t *f, doublewf_t *w ) {

  int i = 0, k = 0;

  if ( ! f || ! w ) {
    bpm_error( "Invalid pointers in apply_filter()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // reset the buffers
  _reset( f->xv );
  _reset( f->yv );

  _reset( f->xv_ac );
  _reset( f->yv_ac );

  if( f->options & FIR ) {
    // apply filter for Finite Impulse response
    
    // pre-fill anti-causal coefficient array
    if ( f->options & ANTICAUSAL ) {
      for ( k=0; k<f->nxc_ac; k++ ) f->xv_ac[k] = w->wf[k] / f->gain;
    }

    for ( i=0; i<f->ns; i++ ) {

      // if the filter is causal :

      if ( f->options & CAUSAL ) {
	_shift_down( f->xv, f->nxc );
	f->xv[f->nxc-1] = w->wf[i] / f->gain;
	f->wfbuffer[i] = 0.;
	for( k=0; k<f->nxc; k++ ) f->wfbuffer[i] += ( f->xc[k] * f->xv[k] ); 
      }

      if ( f->options & ANTICAUSAL ) {
	// if the filter is anti-causal :
	// The anti-causal coefficients are stored in the xc_ac[0]...xc_ac[nxc_ac-1] array..., 
	if (i>0) _shift_down( f->xv_ac, f->nxc_ac );

	f->xv_ac[f->nxc_ac-1] = ( ( i+f->nxc_ac <= f->ns ) ? w->wf[i+f->nxc_ac-1] / f->gain : 0. );

	for ( k=1; k<f->nxc_ac; k++ ) f->wfbuffer[i] += ( f->xc_ac[k] * f->xv_ac[k] );

      }
    }

  } else {
    // apply filter for Infinite Impulse response
    
    for ( i=0; i<f->ns; i++ ) {

      _shift_down( f->xv, f->nxc );
      f->xv[f->nxc-1] = w->wf[i] / f->gain;

      _shift_down( f->yv, f->nyc );
      f->yv[f->nyc-1] = 0.;
      
      // the x contrib
      for ( k=0; k < f->nxc;   k++ ) f->yv[f->nyc-1] += ( f->xc[k] * f->xv[k] );

      // the y contrib
      for ( k=0; k < f->nyc-1; k++ ) f->yv[f->nyc-1] += ( f->yc[k] * f->yv[k] );

      f->wfbuffer[i] = f->yv[f->nyc-1];
    }


  }

  // copy buffer to input waveform
  for(i=0;i<f->ns;i++) w->wf[i] = f->wfbuffer[i];

  return BPM_SUCCESS;
}


