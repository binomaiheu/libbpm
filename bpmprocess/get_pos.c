/**
   @file
   @ingroup processing
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int get_pos( double Q, double I, double IQphase, double posscale,
	     double *pos ) {

  *pos = -DBL_MAX;

  if ( !pos ) {
    bpm_error( "Invalid pointer argument in get_pos(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  *pos = posscale * ( I*cos( IQphase ) + Q*sin( IQphase ) );

  return BPM_SUCCESS;
}
