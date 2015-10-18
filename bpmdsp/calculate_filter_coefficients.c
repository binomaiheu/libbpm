/** 
    @file
    @ingroup dsp
*/

#include "bpm/bpm_dsp.h"

int _expand_complex_polynomial( complex_t *w, int n, complex_t *a ) {
  
  /**
     Calculate the polynomial coefficients in 
        a0 + a1 * z + a2 * z^2 + a3 * z^3 + ...  = (z-w1)(z-w2)(z-w3)...
     from the n polynomial's zero's "w"  
     returns the results in a, the array of coefficients...
  */
  int i = 0, j = 0;

  if ( !w || !a ) {
    bpm_error( "Invalid pointers in _expand_complex_polynomial", __FILE__, __LINE__ );
    return 1;
  }

  a[0] = complex( 1., 0. );
  for ( i=0; i<n; i++ ) a[i+1] = complex( 0., 0. ); 
  for ( i=0; i<n; i++ ) {
    for ( j=n; j>=1; j-- ) a[j] = c_sum( c_mult( c_neg(w[i]), a[j] ), a[j-1]);
    a[0] = c_mult( a[0], c_neg( w[i] ) );
  }

  // check whether a[] are all real numbers 
  for ( i=0; i<n+1; i++ ) {

    if ( fabs( c_imag( a[i]) ) > FILT_EPS ) {
      bpm_error( "Poles/zeros not complex conjugates", __FILE__, __LINE__ );
      return 1;
    }
  }

  return 0;
}

// ----------------------------------------------------------------------------

complex_t _eval_complex_polynomial( complex_t *a, int n, complex_t z ) {

  complex_t p = complex( 0., 0. );
  int i = 0;
  for ( i=n; i>=0; i-- ) p = c_sum( c_mult( p, z ), a[i] );

  return p;
}


// ----------------------------------------------------------------------------

int calculate_filter_coefficients( filter_t *f ) {

  /**
     Calculates the filter coefficients from the poles and zeros in the 
     cplane representation...
     Also calculates the filter gains...
  */
  complex_t topco[MAXPZ+1], botco[MAXPZ+1];
  double theta = 0;
  int i, nz, np;

  nz = f->cplane->nzeros;
  np = f->cplane->npoles;

  if ( _expand_complex_polynomial( f->cplane->zero, nz, topco ) ) return 1;
  f->nxc = nz+1;

  if ( _expand_complex_polynomial( f->cplane->pole, np, botco ) ) return 1;
  f->nyc = np+1;
  
  // calculate gain of filter at DC, HF and FC
  f->dc_gain = c_div( _eval_complex_polynomial( topco, nz, complex( 1., 0. ) ),
		      _eval_complex_polynomial( botco, np, complex( 1., 0. ) ) );

  theta = 2.*PI*0.5*( f->alpha1 + f->alpha2 );
  f->fc_gain = c_div( _eval_complex_polynomial( topco, nz, complex(cos(theta),sin(theta)) ), 
		      _eval_complex_polynomial( botco, np, complex(cos(theta),sin(theta)) ) );

  f->hf_gain = c_div( _eval_complex_polynomial( topco, nz, complex( -1., 0.) ),
		      _eval_complex_polynomial( botco, np, complex( -1., 0.) ) );


  // select the filter gain 
  if ( f->options & LOWPASS )      f->gain = c_abs( f->dc_gain );
  else if( f->options & HIGHPASS ) f->gain = c_abs( f->hf_gain );
  else if( f->options & BANDPASS ) f->gain = c_abs( f->fc_gain );
  else if( f->options & BANDSTOP ) f->gain = c_abs( c_sqrt( c_mult( f->dc_gain, 
								    f->hf_gain ) ) );
  else f->gain = 1.;

  // the coefficients
  for ( i=0; i<f->nxc; i++ ) f->xc[i] =  c_real( topco[i] ) / c_real( botco[f->nyc-1] );
  for ( i=0; i<f->nyc; i++ ) f->yc[i] = -c_real( botco[i] ) / c_real( botco[f->nyc-1] );

  return 0;
}
