/**
   @file
   @ingroup processing
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int get_slope( double Q, double I, double IQphase, double slopescale,
	       double *slope ) {

  *slope = -DBL_MAX;
  if ( !slope ) {
    bpm_error( "Invalid pointer argument in get_slope(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  *slope = slopescale * ( -I*sin(IQphase) + Q*cos(IQphase) );

  return BPM_SUCCESS;
}
