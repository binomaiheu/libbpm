/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Find the median value of the given array. Basically a wrapper for nr_select
   
   @return The value of the median element
*/
double nr_median( int n, double *arr ) {
  
  if( ! arr ) {
    bpm_error( "Invalid array in nr_median(...)", 
	       __FILE__, __LINE__ );
    return -DBL_MAX;
  }

  return nr_select( (int) (n / 2), n, arr );
}
	     
