/**
  @file 
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Sample from a uniform distribution between (and exluding) the 
   upper and lower values.
   
   @param lower the lower range for the generation
   @param upper the upper range for the generation
   
   @return the value of the uniform deviate, returns -DBL_MAX if the seed
           was not set correctly before
*/
double nr_ranuniform( double lower, double upper ) {

  double res = -DBL_MAX;
  
  res = nr_ran1( &bpm_rseed );
  res = res * ( upper - lower ) + lower;
  
  return res;
}
