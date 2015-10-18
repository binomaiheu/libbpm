/**
   @file
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int time_to_sample( double fs, int ns, double t, int *iS ) {
  /*
    @ingroup processing
    
    @param  fs       sampling frequency
    @param  ns       number of samples
    @param  t  the queried time sample
    @param  iS  the returned sample number

    @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  double ti;
  int i;

  if ( t < 0 ) {
    *iS = 0;
    return BPM_SUCCESS;
  }
  
  for ( i=0; i<ns; i++ ) {
    ti = (double) i / fs;
    if ( ABS( t - ti ) <  ( 1./fs ) ) break;
  }

  *iS = i;

  
  return BPM_SUCCESS;
}
