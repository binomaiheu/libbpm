/**
   @file
   @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Returns the incomplete gamma function NR C6.2, p219
*/
int nr_gcf(double *gammcf, double a, double x, double *gln) {

  int i;
  double an,b,c,d,del,h;

  *gln = nr_gammln(a);
  if( *gln == -DBL_MAX ) {
    bpm_error( "nr_gammln failed in nr_gcf(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  b=x+1.0-a;
  c=1.0/GCF_FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=GCF_ITMAX;i++)
    {
      an = -i*(i-a);
      b +=2.0;
      d=an*d+b;
      if (ABS(d) < GCF_FPMIN) d=GCF_FPMIN;
      c=b+an/c;
      if (ABS(c) < GCF_FPMIN) c=GCF_FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;
      if (ABS(del-1.0) < GCF_EPS) break;
    }

  if (i > GCF_ITMAX) {
    bpm_warning( "A too large, GCF_ITMAX too small in nr_gcf(...)", 
		   __FILE__, __LINE__ );
  }

  *gammcf = exp(-x+a*log(x)-(*gln))*h;

  return BPM_SUCCESS;
}
