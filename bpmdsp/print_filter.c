/**
   @file
   @ingroup dsp
*/

#include "bpm/bpm_dsp.h"

void print_filter( FILE *of, filter_t *f ) {

  int i = 0;

  if ( ! of || ! f ) {
    bpm_error( "Invalid pointer in print_filter()", __FILE__, __LINE__  );
    return;
  }

  fprintf( of, "Filter: %s\n", f->name );

  // print some more bla bla here if i find the time...

  if(f->cplane) print_filter_representation( of, f->cplane );
  fprintf( of, "\n" );

  fprintf( of, " - filter gains: %s\n", f->name );
  fprintf( of, "   DC mag= %f, phase= %f pi\n", c_abs(f->dc_gain), c_arg(f->dc_gain)/PI );
  fprintf( of, "   FC mag= %f, phase= %f pi\n", c_abs(f->fc_gain), c_arg(f->fc_gain)/PI );
  fprintf( of, "   HF mag= %f, phase= %f pi\n", c_abs(f->hf_gain), c_arg(f->hf_gain)/PI );
  fprintf( of, "   Filter gain = %f\n", f->gain );
  fprintf( of, "\n" );

  fprintf( of, " - Recurrence relation :\n" );
  fprintf( of, "   y[n] = \n" );

  if ( f->options & CAUSAL ) { // print the causal part
    for ( i=0; i<f->nxc; i++ ) {
      if ( fabs( f->xc[i]) > FILT_EPS ) {
	fprintf( of, "        %s %14.10f * x[n-%d]\n", 
		 f->xc[i]>=0.?"+":"-", fabs(f->xc[i]), f->nxc-1-i );
      }
    }
  }


  if ( f->options & ANTICAUSAL ) { // print the anti-causal part
    for ( i=1; i<f->nxc_ac; i++ ) {
      if ( fabs( f->xc_ac[i]) > FILT_EPS ) {
	fprintf( of, "        %s %14.10f * x[n+%d]\n", 
		 f->xc_ac[i]>=0.?"+":"-", fabs(f->xc_ac[i]), i );
      }
    }
  }


  for ( i=0; i<f->nyc-1; i++ ) {
    if ( fabs(f->yc[i]) > FILT_EPS ) {
      fprintf( of, "        %s %14.10f * y[n-%d]\n", f->yc[i]>=0.?"+":"-", 
               fabs(f->yc[i]), f->nyc-1-i );
    }
  }

  fprintf( of, "\n" );

  return;
}
