/**
   @file
   @ingroup dsp
*/

#include "bpm/bpm_dsp.h"

int filter_step_response( filter_t *f, doublewf_t *w, int itrig ) {

  /**
     Produces a stepresponse for the filter, step is defined by the trigger sample number
     the starting level and the endlevel
  */

  int i = 0;

  if ( ! w || ! f ) {
    bpm_error( "Invalid pointers in filter_step_reponse(...)", __FILE__, __LINE__ );
    return 1;
  } 

  // fill the step
  for ( i = 0; i < f->ns; i++ ) w->wf[i] = i < itrig ? 0. : 1.;

  // apply the filter to the step function
  if ( apply_filter( f, w ) ) {
    bpm_error( "Unable to apply filter in filter_step_response(...)",
	       __FILE__, __LINE__ );
    return 1;
  }

  return 0;
}
