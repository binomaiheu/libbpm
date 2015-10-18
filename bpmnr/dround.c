/**
   @file
   @ingroup nr
*/

double dround( double x ) {
  if( ( x - (int) x ) >= 0.5 ) return (double) ( (int) x + 1. );
  return (double) ( (int) x );
}
