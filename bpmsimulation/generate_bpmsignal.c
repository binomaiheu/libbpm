/**
   @file
   @ingroup sim
*/

#include <bpm/bpm_simulation.h>
#include <bpm/bpm_wf.h>

int generate_bpmsignal( bpmconf_t *bpm, bpmmode_t *mode, beamconf_t *beam, doublewf_t *rf ) {

  int j,k = 0;
  int shift = 0;
  static char msg[80];
  //doublewf_t *rf_bunch = NULL;

  if ( ! bpm || !mode || ! beam || ! rf ) {
    bpm_error( "Invalid pointer arguments in generate_bpmsignal(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! mode->response ) {
    bpm_error( "Mode response is not defined in generate_bpmsignal(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  //printf("Mode response has %i samples\n", mode->response->ns);

  if ( ! mode->buffer ) {
    mode->buffer = doublewf( mode->response->ns, mode->response->fs );
  }  

  //if ( ! rf_bunch ) {
  //  bpm_error( "Can not allocate a temporary array in generate_bpmsignal(...)",
  //              __FILE__, __LINE__ );
  //  return BPM_FAILURE;
  //}

  //printf("Processing %i bunches...\n", beam->nbunches);

  // loop over all bunches
  for ( j=0;j<beam->nbunches;j++ ) {

    //printf("Processing bunch %i\n", j);

    doublewf_reset( mode->buffer );

    if ( add_mode_response( bpm, mode, &( beam->bunch[j] ), mode->buffer ) ) {
      sprintf( msg, "Cannot add response for mode %s in generate_bpmsignal", 
               mode->name );
      bpm_error( msg, __FILE__, __LINE__ );
      return BPM_FAILURE;
    }

    //printf( "Adding signal for bunch %i arriving at dt = %f ns\n", j, beam->bunch[j].arrival_time/_nsec__ );

    // find the first sample of the multibunch waveform
    // affected by the current bunch
    if ( beam->bunch[j].arrival_time == 0. ) {
      shift = 0;
    } else {
      shift = floor( beam->bunch[j].arrival_time * rf->fs ) + 1;
    }

    if ( shift < 0 || shift > rf->ns ) {
      bpm_error( "Sorry, but I've gone mad in generate_bpmsignal(...)",
                __FILE__, __LINE__ );
      return BPM_FAILURE;
    }

    //printf("That means it's shifted by %i samples\n", shift);

    k = 0; // reset the sample counter

    // add the bunch produced waveform to the multibunch waveform
    // interpolating inside the waveform to get the timing right
    while ( shift+k <= rf->ns && k/rf->fs < mode->buffer->ns/mode->buffer->fs ) {

      rf->wf[shift+k] += doublewf_getvalue( mode->buffer,
			 (shift+k)/rf->fs - beam->bunch[j].arrival_time,
			 WF_QUADRATIC );

		         //WF_LANCZOS );

      k++;

    }  /* while ( shift + k <= rf->ns && k <= rf_bunch->ns ) */

    //printf("k = %i\n", k);

    //printf("*");

  } /* for ( j=0;j<beam->nbunches;j++ ) */

  //printf(" Done\n");

  //doublewf_delete( rf_bunch );

  //} /* for ( i=0; i<bpm->nmodes; i++ ) */ 

  return BPM_SUCCESS;
}
