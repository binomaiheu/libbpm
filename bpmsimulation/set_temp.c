/**
   @file
   @ingroup sim
*/
#include <bpm/bpm_interface.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_simulation.h>

double  ambient_temp = 293.;

/**
   Sets up the ambient temperature
   @param TK ambient temperature in Kelvin
   
   @returns BPM_SUCCESS
*/
int set_temp( double TK ) {

  ambient_temp = TK;

  return BPM_SUCCESS;
}
// end of file
