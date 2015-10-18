/**
   @file
   @ingroup dsp
*/

#include "bpm/bpm_dsp.h"

void print_filter_representation( FILE* of, filterrep_t *r ) {
  
  /**
     Display filter representation
  */
  int i = 0;

  if ( ! of || ! r ) return;

  fprintf( of, " - filter zeros : %d \n", r->nzeros );
  for ( i = 0; i< r->nzeros; i++ ) {
    fprintf( of, "   z[%d] = %14.10f %s %14.10f * i\n", i, c_real( r->zero[i]), 
	     ( c_imag( r->zero[i] ) < 0. ) ? "-" : "+", fabs( c_imag( r->zero[i]) ) );
  }

  fprintf( of, " - filter poles : %d \n", r->npoles );
  for ( i = 0; i< r->npoles; i++ ) {
    fprintf( of, "   p[%d] = %14.10f %s %14.10f * i\n", i, c_real( r->pole[i]), 
	     ( c_imag( r->pole[i] ) < 0. ) ? "-" : "+", fabs( c_imag( r->pole[i]) ) );
  }
  
  fflush( of );

  return;
}
