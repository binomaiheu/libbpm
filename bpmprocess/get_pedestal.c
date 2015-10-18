/**
   @file
   @ingroup processing
*/
#include <math.h>

#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int get_pedestal( doublewf_t *wf, int range, double *offset, double *rms ) {

  wfstat_t s;
  int i;

  if ( ! wf || ! offset ) {
    bpm_error( "Invalid pointer argument in get_pedestal(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // get basic statistics from the waveform in the beginning
  if ( doublewf_basic_stats( wf, 0, range, &s ) == BPM_FAILURE ) {
    bpm_error( "Error retreiving basic stats in get_pedestal()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  *offset = s.mean;
  *rms    = s.rms;

  return BPM_SUCCESS;
}
