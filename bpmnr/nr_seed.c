/** 
  @file 
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>


long bpm_rseed = -15041979; /**< the global random seed variable */


/**
   Set the random seed 'idum' to enable other random functions to work
   
   @param seed a random seed
   
   @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
*/
int nr_seed( long seed ) {

  if ( seed == 0 ) {
    bpm_error("Cannot have a 0 random seed", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  bpm_rseed = ( ( seed > 0 ) ? -seed : seed );
  
  return BPM_SUCCESS;
}
