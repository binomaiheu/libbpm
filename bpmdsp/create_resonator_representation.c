/**
   @file
   @ingroup dsp
*/
#include "bpm/bpm_dsp.h"


complex_t _reflect(complex_t z) { 
  double r = c_abs(z);
  return c_div( z, complex( r*r, 0. ) );
}

// -------------------------------------------------------------------------

filterrep_t* create_resonator_representation( filter_t *f ) {

  double       theta, r, thm, th1, th2, phi;
  int          i, cvg;
  complex_t    z, g;
  complex_t    topco[MAXPZ+1], botco[MAXPZ+1];
  filterrep_t *zrep;

  // allocate memory for the filterrepresentation...
  zrep = (filterrep_t*) calloc( 1, sizeof( filterrep_t) );

  if ( ! zrep ) {
    bpm_error( "Cannot allocate memory for resonator representation.", __FILE__, __LINE__ );
    return NULL;
  }

  // for all : 
  zrep->nzeros  = 2;
  zrep->npoles  = 2;
  zrep->zero[0] = complex(  1., 0. );
  zrep->zero[1] = complex( -1., 0. );
  theta         = 2.*PI*f->alpha1;

  if ( f->Q <= 0. ) {
    // negative Q, so assume infinte Q : pure oscillator, calculate the poles
    z = c_exp( complex(0., theta) );
    zrep->pole[0] = z;
    zrep->pole[1] = c_conj( z );
  } else {
    // finite Q factor for the resonator, calculate the poles
    _expand_complex_polynomial( zrep->zero, zrep->nzeros, topco);
    r = exp(-theta / (2.0 * f->Q));
    thm = theta;
    th1 = 0.0;
    th2 = PI;
    cvg = 0;
    for ( i=0; ( i<MAX_RESONATOR_ITER ) && ( cvg == 0 ); i++ ) {
      z = complex( r*cos(thm), r*sin(thm));
      zrep->pole[0] = z; 
      zrep->pole[1] = c_conj(z);

      _expand_complex_polynomial( zrep->pole, zrep->npoles, botco);

      g = c_div( _eval_complex_polynomial( topco, zrep->nzeros, complex(cos(theta),sin(theta))), 
		 _eval_complex_polynomial( botco, zrep->npoles, complex(cos(theta),sin(theta))) );
      phi = g.im / g.re;
      if ( phi > 0.0 ) th2 = thm; else th1 = thm;
      if ( fabs(phi) < FILT_EPS ) cvg = 1;
      thm = 0.5 * ( th1 + th2 );
    }
    if ( ! cvg ) {
      bpm_error( "Finite Q resonator failed to converge on pole/zero calculation.",
		 __FILE__, __LINE__ );
      free(zrep);
      return NULL;
    }

  }


  // adjust the zeros for a bandstop resonator
  if ( f->options & BANDSTOP ) {
    theta = 2.*PI*f->alpha1;
    z = complex( cos(theta), sin(theta) );
    zrep->zero[0] = z;
    zrep->zero[1] = c_conj( z );
  }

  // adjust the zeros for an allpas resonator
  if ( f->options & ALLPASS ) {
    zrep->zero[0] = _reflect( zrep->pole[0] );
    zrep->zero[1] = _reflect( zrep->pole[1] );
  }


  return zrep;
}
