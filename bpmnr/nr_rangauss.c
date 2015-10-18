/**
  @file 
  @ingroup nr
*/
#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Sample a given Gaussian distribution using ran1 as the source of the
   uniform deviate between 0 and 1. Nicked from numerical recipes, C7.2, p289
   
   @param mean the mean of the gaussian
   @param std_dev the standard deviation of the gaussian
   
   @return a gaussian deviate, returns -DBL_MAX if the random seed is not
   set properly before
*/
double nr_rangauss( double mean, double std_dev ) {

  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  double val = -DBL_MAX;

  if ( iset == 0 ) {
    do {
      v1 = 2.0*nr_ran1( &bpm_rseed ) - 1.0;
      v2 = 2.0*nr_ran1( &bpm_rseed ) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0 );

    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    val = (v2*fac * std_dev) + mean;
  } else  {
    iset = 0;
    val = (gset * std_dev) + mean;
  }

  return val;
}
  
