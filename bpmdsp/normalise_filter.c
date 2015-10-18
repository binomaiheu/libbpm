/**
   @file
   @ingroup dsp
*/
#include "bpm/bpm_dsp.h"

int normalise_filter( filter_t *f, filterrep_t *s ) {

  double w0, bw;
  complex_t hba, tmp;
  double w1 = 2. * PI * f->w_alpha1;
  double w2 = 2. * PI * f->w_alpha2;
  int i = 0;

  if ( ! f || ! s ) {
    bpm_error( "Invalid pointer in normalise_filter()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // some calculations for bandpass/stop filters
  w0 = sqrt( w1 * w2 );
  bw = w2 - w1;

  if ( f->options & LOWPASS ) {
    // normalise for low pass filter

    for ( i=0; i< s->npoles; i++ ) s->pole[i] = c_scale( w1, s->pole[i] );
    s->nzeros = 0;

  } else if ( f->options & HIGHPASS ) {
    // normalise for high pass filter

    for ( i=0; i<s->npoles; i++ ) s->pole[i] = c_div( complex(w1,0.), s->pole[i] );
    for ( i=0; i<s->npoles; i++ ) s->zero[i] = complex( 0., 0.);
    s->nzeros = s->npoles;

  } else if ( f->options & BANDPASS ) {
    // normalise for band pass filter

    for( i=0; i<s->npoles; i++ ) {
      hba = c_scale( .5 * bw, s->pole[i] );
      
      tmp = c_sqrt( c_diff( complex( 1., 0.), 
			    c_div( complex( w0*w0, 0. ), c_sqr( hba ) ) ) );
      
      s->pole[i]           = c_mult( hba,  c_sum( complex( 1., 0.), tmp ) );
      s->pole[s->npoles+i] = c_mult( hba, c_diff( complex( 1., 0.), tmp ) );
    }

    for( i=0; i<s->npoles; i++ ) s->zero[i] = complex( 0., 0.);
    s->nzeros  = s->npoles;
    s->npoles *= 2.;

  } else if ( f->options & BANDSTOP ) {
    // normalise for band stop filter
    for( i=0; i<s->npoles; i++ ) {
      hba = c_div( complex( .5 * bw, 0.), s->pole[i] );
      tmp = c_sqrt( c_diff( complex( 1., 0.), 
			    c_div( complex( w0*w0, 0. ), c_sqr( hba ) ) ) );

      s->pole[i]             = c_mult( hba,  c_sum( complex(1.,0.), tmp ) );
      s->pole[s->npoles+i] = c_mult( hba, c_diff( complex(1.,0.), tmp ) );
    }

    for( i=0; i<s->npoles; i++ ) {
      s->zero[i]             = complex( 0.,  w0 );
      s->zero[s->npoles+i] = complex( 0., -w0 ); 
    }
    s->npoles *= 2.;
    s->nzeros  = s->npoles;

  } 

  return BPM_SUCCESS;
}
