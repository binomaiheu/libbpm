/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>


/**
   Checks whether the given double is an integer value, handy for 
   doing domain checking to prevent e.g. the function nr_gammln print
   out "nan" or "inf" values...
   
   For double precision, this check is accurate to 1.0E-323 ... should be
   enough ;-)
   
   @param x floating point argument
   
   @return TRUE if argument is indeed an integer value, FALSE if not
*/ 
int nr_is_int( double x ) {

  if( ABS( x - (int) x ) == 0. ) return TRUE;
  
  return FALSE;
}

// -----------------------------------------------------------------------------

/**
   Checks whether the input argument is an integer power of 2, like
   256, 1024 etc...
   
   @param n given unsigend long argument for which to check this
   
   @return FALSE if not a power of 2. The routine returns the precise power
   ( > 1 ) if the integer is indeed a power of 2
*/
int nr_is_pow2( unsigned long n ) {

  int p = 0;
  int r = 0;

   do {
     r = n % 2;  n /= 2;  p++;
   } while ( r == 0 && n > 1 );

   if( r != 0 ) return FALSE;

   return p;
}

