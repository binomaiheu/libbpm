/**
   @file
   @ingroup rf
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_rf.h>

int     rf_nsamples   = 65536; 
double  rf_samplefreq = 20. * _GHz__;

/**
   Sets up the sampling of internal RF waveform representation
   @param nsamples the number of samples
   @param sfreq the internal sampling frequency
   
   @returns BPM_SUCCESS
*/
int rf_setup( int nsamples, double sfreq ) {

  
  rf_nsamples   = nsamples;
  rf_samplefreq = sfreq;

  return BPM_SUCCESS;
}
// end of file
