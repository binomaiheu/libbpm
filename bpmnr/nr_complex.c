/**
   @file
   @ingroup nr
*/

#include "bpm/bpm_nr.h"


/* Defines a complex number from the real and imag parts */
complex_t complex( double re, double im ) {
  complex_t z;
  z.re = re;
  z.im = im;
  return z;
}

double c_real( complex_t z ) { return z.re; }
double c_imag( complex_t z ) { return z.im; }
double c_abs( complex_t z ) { return sqrt( z.re*z.re+z.im*z.im ); }
double c_abs2( complex_t z ) { return z.re*z.re+z.im*z.im; }
double c_arg( complex_t z ) { return atan2( z.im, z.re ); }

complex_t c_conj( complex_t z ) { return complex( z.re, -z.im ); }
complex_t c_neg( complex_t z ) { return complex( -z.re, -z.im ); }

complex_t c_sum( complex_t z1, complex_t z2 ) { return complex( z1.re + z2.re, z1.im + z2.im ); }
complex_t c_diff( complex_t z1, complex_t z2 ) { return complex( z1.re - z2.re, z1.im - z2.im ); }

complex_t c_mult( complex_t z1, complex_t z2 ) {
  return complex( z1.re*z2.re - z1.im*z2.im, 
		  z1.re*z2.im + z1.im*z2.re );
}

complex_t c_scale( double r, complex_t z ) { return complex( r*z.re, r*z.im); }


complex_t c_div( complex_t z1, complex_t z2 ) {
  double r = c_norm2( z2 );
  return complex( ((z1.re*z2.re) + (z1.im*z2.im)) / r, 
		  ((z1.im*z2.re) - (z1.re*z2.im)) / r );
}

complex_t c_sqr( complex_t z ) { return c_mult( z, z ); }

double c_norm2( complex_t z ) { return ( z.re*z.re + z.im*z.im ); }

complex_t c_exp( complex_t z ) { return complex( exp(z.re)*cos(z.im), exp(z.re)*sin(z.im) ); }

complex_t c_sqrt( complex_t z ) {
  double    r = c_abs( z );
  complex_t u = complex( sqrt( (r + z.re ) / 2. ), sqrt( ( r - z.re ) / 2. ) );
  if( z.im < 0.0 ) u.im = -u.im;
  return u;
}

/* Is z1 equal to z2 ? */
int c_isequal( complex_t z1, complex_t z2 ) {
  if ( ( z1.re == z2.re ) && ( z1.im == z2.im ) ) {
    return 1;
  } else return 0;
}

