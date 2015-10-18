/**
   @file
   @ingroup dsp
*/

#include "bpm/bpm_dsp.h"

filterrep_t* zplane_transform( filter_t *f, filterrep_t *s ) {

  int i = 0;
  filterrep_t *r;
  
  if ( ! s ) {
    bpm_error( "Invalid pointer argument in zplane_transform(...).", __FILE__, __LINE__ );
    return NULL;
  }


  // allocate memory for the filterrepresentation...
  r = (filterrep_t*) calloc( 1, sizeof( filterrep_t) );

  if ( ! r ) {
    bpm_error( "Cannot allocate memory for z-plane representation.", __FILE__, __LINE__ );
    return NULL;
  }

  // the number of poles and zeros
  r->npoles = s->npoles;
  r->nzeros = s->nzeros;

  if ( f->options & MATCHED_Z_TRANSFORM ) {

    /* Convert s plane poles to z plane using matched z transform */
    for ( i=0; i<r->npoles; i++ ) r->pole[i] = c_exp( s->pole[i] );
    for ( i=0; i<r->nzeros; i++ ) r->zero[i] = c_exp( s->zero[i] );
    
  } else {
    /* Default case... */
    /* Convert s plane poles to z plane using bilinear transform */
    for ( i=0; i<r->npoles; i++ ) { 
      r->pole[i] = c_div( c_sum( complex(2.,0.), s->pole[i] ),
			  c_sum( complex(2.,0.), c_scale(-1.0,s->pole[i])) );
    }
    
    for ( i=0; i<r->nzeros; i++ ) {
      r->zero[i] = c_div( c_sum( complex(2.,0.), s->zero[i] ),
			  c_sum( complex(2.,0.), c_scale(-1.0,s->zero[i])) );
    }
    
    while ( r->nzeros < r->npoles ) r->zero[r->nzeros++] = complex( -1.0, 0. );

  } /* switch ( mode ) */


  return r;
}


