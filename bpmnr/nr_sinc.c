/**
   @file
   @ingroup nr
*/

#include "bpm_nr.h"

double sinc( double x ) {
  if ( fabs(x) < 1.e-11 ) return 1.0;
  return sin(PI*x)/(PI*x);
}

double lanczos( double x, int a ) {
  if ( x ==  0.0 ) return 1.0;
  if ( x <= (double) -a || x >= (double) a ) return 0.0;
  return sinc(x) * sinc( x / (double) a );
}
