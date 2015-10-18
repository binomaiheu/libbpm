/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>


/**
   Replaces data[1..2*nn] by its discrete Fourier transform,
   if isign is input as 1, or replaces data[1..2*nn] by nn times its
   inverse discrete Fourier transform if isign is input as -1.
   
   data is a complex arry of length nn, or equivalently a real array of
   length 2*nn. nn MUST !!! be an integer power of 2, this is not checked
   for...
   
   BM. 15.08.2005... added this check ;-))
   
   Perform an FFT, NR S12.2 pg507
   See: http://www.library.cornell.edu/nr/cbookcpdf.html
   
   
   @param data array with data
   @param nn number of data points, note that the array length has to be
   at least twice this number
   @param isign sign of transform
   
   @return BPM_FAILURE upon failure, BPM_SUCCESS upon success
*/

int nr_four1( double data[], unsigned long nn, int isign ) {

  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  if ( ! nr_is_pow2( nn ) ) {
    bpm_error( "Data length is not power of 2 in nr_four1(...)",
		 __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  n=nn<<1;
  j=1;
  for (i=1;i<n;i+=2)
    {
      if (j > i)
	{
	  SWAP(data[j],data[i]);
	  SWAP(data[j+1],data[i+1]);
	}
      
      m=nn;
      while (m >= 2 && j > m)
	{
	  j -= m;
	  m >>= 1;
	}
      j +=m;
    }

  mmax=2;
  while(n>mmax)
    {
      istep=mmax << 1;
      theta=isign*(2*PI/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      
      for (m=1;m<mmax;m+=2)
	{
	  for (i=m;i<=n;i+=istep)
	    {
	      j=i+mmax;
	      tempr=wr*data[j]-wi*data[j+1];
	      tempi=wr*data[j+1]+wi*data[j];
	      data[j]=data[i]-tempr;
	      data[j+1]=data[i+1]-tempi;
	      data[i] += tempr;
	      data[i+1] += tempi;
	    }
	  wr=(wtemp=wr)*wpr-wi*wpi+wr;
	  wi=wi*wpr+wtemp*wpi+wi;
	}
      mmax=istep;
    }

  return BPM_SUCCESS;
}
