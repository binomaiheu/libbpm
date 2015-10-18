/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_nr.h>

/**
   Random number generator as nicked from numerical recipes, c7.1, p280
   
   @param idum random seed, note that the global seed is set by bpm_rseed
   @return random number between 0 and 1
*/
double nr_ran1( long *idum ) {
  
  int j;
  long k;
  static long iy = 0;
  static long iv[RAN1_NTAB];
  double temp;

  if (*idum <= 0 || !iy)
    {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=RAN1_NTAB+7;j>=0;j--)
	{
	  k=(*idum)/RAN1_IQ;
	  *idum=RAN1_IA*(*idum-k*RAN1_IQ)-RAN1_IR*k;
	  if (*idum < 0) *idum += RAN1_IM;
	  if (j < RAN1_NTAB) iv[j] = *idum;
	}
      iy=iv[0];
    }
  k=(*idum)/RAN1_IQ;
  *idum=RAN1_IA*(*idum-k*RAN1_IQ)-RAN1_IR*k;
  if (*idum < 0) *idum += RAN1_IM;
  j=iy/RAN1_NDIV;
  iy=iv[j];
  iv[j]=*idum;
  if((temp=RAN1_AM*iy) > RAN1_RNMX) return RAN1_RNMX;
  else return temp;
}
