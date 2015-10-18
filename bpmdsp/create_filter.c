/**
   @file
   @ingroup dsp
*/
#include <string.h>
#include <stdlib.h>

#include "bpm/bpm_dsp.h"

filter_t* create_filter( char name[], unsigned int options, int order, int ns,
			 double fs, double f1, double f2, double par ) {
  
  filterrep_t *s;
  int n = 0;

  filter_t *f = (filter_t*) calloc( 1, sizeof(filter_t) );
  if ( ! f ) {
    bpm_error( "Couldn't reserve memory for filter", __FILE__, __LINE__ );
    return NULL;
  }

  strncpy( f->name, name, 79 );
  f->options = options;
  f->order   = order;

  f->ns      = ns;
  f->fs      = fs;
  f->f1      = f1;
  f->f2      = f2;

  f->cheb_ripple  = 0.;
  f->gauss_cutoff = 0.001;
  f->Q            = -1.;  // initialise to infinite Q, so pure oscillator for resonator filter

  if ( f->options & CHEBYSHEV ) {
    if ( par <  0. ) {
      f->cheb_ripple = par;
    } else {
      bpm_warning( "Invalid Chebyshev ripple, setting default !", __FILE__, __LINE__ );
    }
  }

  if ( f->options & RESONATOR ) {
    if ( par > 0. ) {
      f->Q = par;
    } else {
      bpm_warning( "Q factor <= 0, assuming pure oscillator !", __FILE__, __LINE__ );
    }
  }

  if ( f->options & GAUSSIAN ) {
    if ( par <=  0. ) {
      bpm_warning( "Invalid gaussian cutoff, setting default !", __FILE__, __LINE__ );
    } else {
      f->gauss_cutoff = par;
    }
  }

  if ( f->fs > 0. ) {
    f->alpha1   = f->f1 / f->fs;
    f->alpha2   = f->f2 / f->fs;

    if ( f->options & NO_PREWARP ) {
      f->w_alpha1 = f->alpha1;  
      f->w_alpha2 = f->alpha2;  
    } else {
      f->w_alpha1 = tan( PI * f->alpha1 ) / PI;  // pre-warped alpha1
      f->w_alpha2 = tan( PI * f->alpha2 ) / PI;  // pre-warped alpha2
    }
  } else {
    bpm_error( "Invalid sampling frequency in create_filter(...)", __FILE__, __LINE__ );
    free(f);
    return NULL;
  }


  // setup the waveform buffer


  f->wfbuffer = (double*) calloc( ns, sizeof(double)  );
  if ( ! f->wfbuffer ) {
    bpm_error( "Cannot allocate memory for waveform buffer in create_filter()",
	       __FILE__, __LINE__ );
    free(f);
    return NULL;
  }
 

  // calculate poles and zeros where applicable
  

  if ( f->options & ( BUTTERWORTH | BESSEL | CHEBYSHEV ) ) {
    // this filter is causal
    f->options |= CAUSAL;

    // calculate s plane represenation and transform to z plane for
    // butterworth, cheybyshev and bessel filters
    s = create_splane_representation( f );
    normalise_filter( f, s );
    f->cplane = zplane_transform( f, s );
    free(s);
  } 

  if ( f->options & RESONATOR ) {
    // this filter is causal
    f->options |= CAUSAL;

    // set some values for the resonator
    f->alpha2   = f->alpha1;  // needed when we calculate fc_gain
    f->w_alpha2 = f->w_alpha1;

    // Directly create z pole representation for the resonator case
    f->cplane = create_resonator_representation( f );
  }


  // calculate the filter coefficients

  
  if ( f->options & GAUSSIAN ) {
    // this filter is non-causal
    f->options |= NONCAUSAL;

    f->cplane = NULL; // set this to 0, no representation for gaussian filter implementation
    if ( gaussian_filter_coeffs( f ) == BPM_FAILURE ) {
      bpm_error( "Failed to calculate gaussian coefficients",__FILE__, __LINE__ );
      free(f->wfbuffer);
      free(f);
      return NULL;
    };

    // filter implementation is FIR
    f->options |=  FIR;
    f->options &= ~IIR;
  } else {
    // calculate general filter coefficients from the pole/zero representation
    calculate_filter_coefficients( f );

    // set whether filter is FIR/IIR
    while (n < f->cplane->npoles && f->yc[n] == 0.0) n++;
    if ( n >= f->cplane->npoles ) { 
      f->options |=  FIR;
      f->options &= ~IIR;
    } else {
      f->options |=  IIR; 
      f->options &= ~FIR; 
    }
  }
  


  return f;
}
