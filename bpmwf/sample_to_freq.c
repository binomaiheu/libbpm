/**
   @file
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

int sample_to_freq( double fs, int ns, int iS, double *f ) {
  /*
    @ingroup processing

    This routine returns the frequency corresponding to the sample number, note
    that this routine is not aware of the nyquist bands, and just keeps on 
    counting from 0 -> fs. 

    @param  fs  sampling frequency
    @param  ns  number of samples
    @param  iS  the queried sample to get the frequency of
    @param  f   the returned frequency

    @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */

  if ( iS < 0. ) {
    *f = 0.;
    return BPM_SUCCESS;
  }
  
  if ( iS >= ns ) {
    *f = fs; // return sampling freq
    return BPM_SUCCESS;
  }

  *f = (double) iS / ((double)ns) * fs;
  
  return BPM_SUCCESS;
}
