/**
   @file
   @ingroup sim
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_simulation.h>

double system_time = 0.;

/**
   Sets up the system clock
   @param ts current time in seconds
   
   @returns BPM_SUCCESS
*/
int set_time( double ts ) {

  system_time = ts;

  return BPM_SUCCESS;
}
// end of file
