/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Calculates the logaritm of the gamma function
   ln[gamma(xx)]. NR C6.1, p 214
   supposed to be correct to double precision
   
   @param xx the argument
   @return the value of ln[gamma(xx)]
*/
double nr_gammln( double xx ) {

  double x,y,tmp,ser;
  int j;
  static double cof[6]={76.18009172947146, -86.50532032941677,
			24.01409824083091, -1.2317395572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};

  if( xx == 0. ) {
    bpm_error( "Argument is 0 in nr_gammln(...)", __FILE__, __LINE__ );
    return -DBL_MAX;
  }

  // Function also undefined for negative integers
  if( xx < 0. && nr_is_int( xx ) ) {
    bpm_error( "Function domain error for nr_gammln(...)", 
	       __FILE__, __LINE__ );
    return -DBL_MAX;
  }

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  for (j=0; j<=5;j++) ser += cof[j]/++y;

  return -tmp+log(2.5066282746310005*ser/x);
}
