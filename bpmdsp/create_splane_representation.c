/**
   @file
   @ingroup dsp
*/
#include "bpm/bpm_dsp.h"

static complex_t bessel_pole[] = { 
  { -1.00000000000e+00, 0.00000000000e+00}, { -1.10160133059e+00, 6.36009824757e-01},
  { -1.32267579991e+00, 0.00000000000e+00}, { -1.04740916101e+00, 9.99264436281e-01},
  { -1.37006783055e+00, 4.10249717494e-01}, { -9.95208764350e-01, 1.25710573945e+00},
  { -1.50231627145e+00, 0.00000000000e+00}, { -1.38087732586e+00, 7.17909587627e-01},
  { -9.57676548563e-01, 1.47112432073e+00}, { -1.57149040362e+00, 3.20896374221e-01},
  { -1.38185809760e+00, 9.71471890712e-01}, { -9.30656522947e-01, 1.66186326894e+00},
  { -1.68436817927e+00, 0.00000000000e+00}, { -1.61203876622e+00, 5.89244506931e-01},
  { -1.37890321680e+00, 1.19156677780e+00}, { -9.09867780623e-01, 1.83645135304e+00},
  { -1.75740840040e+00, 2.72867575103e-01}, { -1.63693941813e+00, 8.22795625139e-01},
  { -1.37384121764e+00, 1.38835657588e+00}, { -8.92869718847e-01, 1.99832584364e+00},
  { -1.85660050123e+00, 0.00000000000e+00}, { -1.80717053496e+00, 5.12383730575e-01},
  { -1.65239648458e+00, 1.03138956698e+00}, { -1.36758830979e+00, 1.56773371224e+00},
  { -8.78399276161e-01, 2.14980052431e+00}, { -1.92761969145e+00, 2.41623471082e-01},
  { -1.84219624443e+00, 7.27257597722e-01}, { -1.66181024140e+00, 1.22110021857e+00},
  { -1.36069227838e+00, 1.73350574267e+00}, { -8.65756901707e-01, 2.29260483098e+00} };


void _add_splane_pole( filterrep_t *r, complex_t z ) {
  if ( c_real( z ) < 0.0 ) r->pole[r->npoles++] = z;
  return;
}

// ---------------------------------------------------------------------------------------

filterrep_t* create_splane_representation( filter_t *f ) {

  char msg[80];
  int p = 0, i = 0;
  double theta, rip, eps, y;
  filterrep_t *r;

  // allocate memory for the filterrepresentation...
  r = (filterrep_t*) calloc( 1, sizeof( filterrep_t) );

  if ( ! r ) {
    bpm_error( "Cannot allocate memory for s-plane representation.", __FILE__, __LINE__ );
    return NULL;
  }

  // initialise number of poles to 0
  r->npoles = 0;

  /* s plane poles for Bessel filter */
  if ( f->options & BESSEL ) {
    p = ( f->order * f->order ) / 4;
    if ( f->order & 1 ) _add_splane_pole( r, bessel_pole[p++] );
    for ( i=0; i<f->order/2; i++ ) {
      _add_splane_pole( r,  bessel_pole[p]        );
      _add_splane_pole( r,  c_conj(bessel_pole[p])  );
      p++;
    }
  }

  /* s plane poles for Chebyshev or Butterworth filter */
  if ( f->options & ( BUTTERWORTH | CHEBYSHEV ) ) {
    for( i=0; i<2*f->order; i++ ) {
      theta = (f->order & 1 ) ? i*PI/f->order : (i+.5)*PI/f->order;
      _add_splane_pole( r, c_exp(complex(0.,theta)) );
    }
  }

  /* Modify for Chebyshev filter */
  if ( f->options & CHEBYSHEV ) {
    if ( f->cheb_ripple >= 0.0 ) {
      bpm_error(  "Chebyshev ripple is must be < 0 dB!", __FILE__, __LINE__ );
      return NULL;
    }

    rip = pow( 10.0, -f->cheb_ripple / 10.0 );
    eps = sqrt( rip - 1.0 );
    y   = asinh( 1.0 / eps ) / (double) f->order;
    if ( y <= 0. ) {
      sprintf( msg, "Chebyshev ripple coefficient is %f, must be > 0", y );
      bpm_error( msg, __FILE__, __LINE__ );
      return NULL;
    } else {
      // modify the exising circular oriented poles...
      for( i=0; i<r->npoles; i++ ) {
        r->pole[i] = complex( c_real(r->pole[i])*sinh(y), c_imag(r->pole[i])*cosh(y) );
      }
    }
  }

  return r;
}
