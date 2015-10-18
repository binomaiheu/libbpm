/**
   @file
   @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Returns the incomplete gamma function. 
   From numerical recipes, C6.2, p218
   
   @return -DBL_MAX upon failure 
*/
double nr_gammq( double a, double x ) {

  double gamser, gammcf, gln;
  double res = -DBL_MAX;

  if ( x < 0.0 || a <= 0.0 ) {
    bpm_error( "Invalid arguments in nr_gammq(...)",
	       __FILE__, __LINE__ );
    return -DBL_MAX;
  }

  if (x < (a+1.0) ) {

    nr_gser( &gamser, a, x, &gln );
    res =  1.0 - gamser;
      
  } else {
    nr_gcf( &gammcf, a, x, &gln );
    res = gammcf;
  }

  return res;
}
