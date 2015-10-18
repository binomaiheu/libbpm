/**
   @file
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int sample_to_time( double fs, int ns, int iS, double *t ) {
  /*
    @ingroup processing
    
    @param  fs       sampling frequency
    @param  ns       number of samples
    @param  t  the queried sample
    @param  iS  the returned sample time

    @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */

  if ( iS < 0. ) {
    *t = 0.;
    return BPM_SUCCESS;
  }
  
  if ( iS >= ns ) {
    *t = (double) ns / fs; // return time of last sample + 1
    return BPM_SUCCESS;
  }

  *t = (double) iS / fs;
  
  return BPM_SUCCESS;
}
