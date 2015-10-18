/**
   @file
   @ingroup sim
*/

#include <bpm/bpm_simulation.h>
#include <bpm/bpm_nr.h>
#include <bpm/bpm_interface.h>

int add_mode_response( bpmconf_t *bpm, bpmmode_t *mode, bunchconf_t *bunch, doublewf_t *rf ) {

  int i,imax = 0;
  complex_t amp;

  if ( ! rf  ) {
    bpm_error( "BPM signal waveform is not pre-allocated in add_mode_response()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // now calculate the complex_t scale factor from the sensitivities
  amp = get_mode_amplitude( bpm, mode, bunch );

  // set the loop length
  if ( rf->ns < mode->response->ns ) imax = rf->ns; else imax = mode->response->ns;

  // scale the response and add it to the signal waveform
  // for efficiency don't use the library functions
  // for the dipole mode add the incline+tilt signal as well
  if (mode->order == 1) {
    for (i=0;i<imax;i++){
      rf->wf[i] += mode->response->wf[i].re*amp.re + mode->response->wf[i].im*amp.im;
    }
  } else {
    for (i=0;i<imax;i++){
      rf->wf[i] += mode->response->wf[i].re*amp.re;
    }
  }
  return BPM_SUCCESS;
}
