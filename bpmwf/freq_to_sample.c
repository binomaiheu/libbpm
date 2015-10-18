/**
   @file
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int freq_to_sample( double fs, int ns, double f, int *iS ) {
  /*
    @ingroup processing
    
    This routine returns the sample number corresponding to the frequency, note
    that this routine is not aware of the nyquist bands, and just keeps on 
    counting from 0 -> fs. 

    @param  fs       sampling frequency
    @param  ns       number of samples
    @param  tsample  the queried frequency of the sample
    @param  isample  the returned sample number

    @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  double fi;
  int i;

  if ( f < 0. ) {
    *iS = 0;
    return BPM_SUCCESS;
  }
  
  for ( i=0; i<ns; i++ ) {
    // get the frequency corrsponding to sample i
    fi = (double) i / ( (double) ns ) * fs;
    if ( ABS( f - fi ) <  ( fs / (double) ns ) ) break;
  }

  *iS = i;

  
  return BPM_SUCCESS;
}
