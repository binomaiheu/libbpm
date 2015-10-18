/**
   @file
   @ingroup processing
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int get_IQ( double amp, double phase, double refamp, double refphase,
	    double *Q, double *I ) {

  if( ! Q || ! I ) {
    bpm_error( "Invalid pointer arguments in get_IQ(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  *I = -DBL_MAX;
  *Q = -DBL_MAX;

  if( refamp > 0. ) {
    *I = amp / refamp * cos( phase - refphase );
    *Q = amp / refamp * sin( phase - refphase );
  } else {
    bpm_warning( "Reference amplitude is 0 in get_IQ(...)",
		 __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}
