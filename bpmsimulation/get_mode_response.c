/**
   @file
   @ingroup sim
*/

#include <bpm/bpm_simulation.h>
#include <bpm/bpm_wf.h>
#include <bpm/bpm_dsp.h>
#include <bpm/bpm_interface.h>

int get_mode_response( bpmmode_t *mode) {

  filter_t *resonator  = NULL;
  doublewf_t *re       = NULL;

  if ( ! mode->response  ) {
    bpm_error( "Buffer for storing the mode response is not defined in add_mode_response()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  complexwf_reset( mode->response );

  re = doublewf( mode->response->ns, mode->response->fs );

  if ( ! re  ) {
    bpm_error( "Cannot allocate memory for a doublewf in add_mode_response()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // define excitation compensating for the energy wasted
  // by filtering in the next step to get amp = 1
  re->wf[0] = mode->Q*mode->response->fs/(2*PI*mode->frequency); 

  //printf("fs =  %f, ns = %i, f0 = %f, Q = %f\n",
  //       mode->response->fs, mode->response->ns,
  //       mode->frequency, mode->Q);

  // mode resonator
  resonator = create_filter( "resonator", RESONATOR | BANDPASS, 0,
			     mode->response->ns, mode->response->fs,
                             mode->frequency, 0., mode->Q );
  
  // apply the filter to get the time domain response
  apply_filter( resonator, re );

  // put the generated response into the buffer
  complexwf_setreal( mode->response, re );

  if ( mode->order == 1 ) {

    //integrate to get a 90 degree phase shift
    doublewf_integrate( re );
    doublewf_scale( 2*PI*mode->frequency, re );

    complexwf_setimag( mode->response, re );

  }

  delete_filter( resonator );
  doublewf_delete( re );

  return BPM_SUCCESS;
}
