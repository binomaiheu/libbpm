/**
   @file
   @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Returns incomplete gamma function. NR 6.2, 218
*/
int nr_gser(double *gamser, double a, double x, double *gln ) {

  int n;
  double sum, del, ap;
  
  // inhibit division by 0 later on
  if( a == 0. ) {
    bpm_error( "a equals 0 in nr_gser(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  *gln = nr_gammln(a);
  if( *gln = -DBL_MAX ) {
    bpm_error( "nr_gammln failed in nr_gser(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  if (x <= 0.0)  {
    if (x < 0.0 ) 
      bpm_warning("x less than 0 in routine nr_gser(...)",  __FILE__, __LINE__ );
    
    *gamser = 0.0;
    return BPM_SUCCESS;
  } else  {
    ap = a;
    del = sum = 1.0 / a;
    for (n=1;n<=GSER_ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (ABS(del) < ABS(sum) * GSER_EPS ) {
	// converged
	*gamser = sum*exp(-x+a*log(x)-(*gln));
	return BPM_SUCCESS;
      }
    } /* for (n=1;n<=GSER_ITMAX;n++) */

    bpm_error( "a too large, GSER_ITMAX too small in nr_gser(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // should not get here
  return BPM_FAILURE;
}
