/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>


/**
   Fit data to a straight line. Nicked from numerical recipes, C15.2, p665
   See: http://www.library.cornell.edu/nr/cbookcpdf.html
   
   @param x array with x values
   @param y array with corresponding y values
   @param ndata number of datapoints
   @param sig array with errors on y datapoints
   @param mwt used weighted (so including errors on datapoints ?)
   @param a fitted slope
   @param b fitted intercept
   @param siga error on fitted slope
   @param sigb error on fitted intercept
   @param chi2 chi2 of fit
   @param q quality factor of fit
   
   @return BPM_FAILURE upon failure, BPM_SUCCESS upon success
*/
int nr_fit ( double *x, double y[], int ndata, double sig[], int mwt, 
	     double *a, double *b, double *siga, double *sigb, 
	     double *chi2, double *q ) {
  int i;
  double wt, t, sxoss, sx=0.0, sy=0.0, st2=0.0, ss, sigdat;

  if( ! ( x && y ) ) {
    bpm_error( "Invalid arguments in nr_fit(...)", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if( ( mwt ) && ! sig ) {
    bpm_error( "Invalid arguments using sig[] in nr_fit(...)", 
	    __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if( ndata < 3 ) {
    bpm_error( "Number of datapoints to small (<3) in nr_fit(...)", 
	    __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  *b=0.0;
  if ( mwt ) {
    ss=0.0;

    for ( i=0; i<ndata;i++ ) {
      if( ABS( sig[i] ) > 0. ) {
	wt=1.0/SQR(sig[i]);
      } else {
	bpm_error( "sig[] contains 0 values in nr_fit(...)",
		     __FILE__, __LINE__ );
	return BPM_FAILURE;
      }
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }

  } else {

    for (i=0; i<ndata;i++) {
      sx += x[i];
      sy += y[i];
    }

    ss=ndata;
  }  /* if ( mwt ) */
  
  if( ABS(ss) > 0. ) {
    sxoss=sx/ss;
  } else {
    bpm_error( "ss is zero in nr_fit(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( mwt ) {

    for (i=0; i<ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      *b += t*y[i]/sig[i];
    }

  } else {
    for (i=0; i<ndata; i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  } /* if ( mwt ) */

  if( ABS(st2) > 0. ) {
    *b /= st2;
    *a = (sy-sx*(*b))/ss;
    *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    *sigb=sqrt(1.0/st2);
  } else {
    bpm_error( "st2 is zero in nr_fit(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  *chi2=0.0;
  *q=1.0;

  if ( mwt == 0 ) {
    for (i=0;i<ndata;i++)
      *chi2 +=SQR(y[i]-(*a)-(*b)*x[i]);

    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;

  } else {

    for (i=0; i<ndata;i++) *chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
    
    *q = nr_gammq(0.5*(ndata-2), 0.5*(*chi2));
    
  } /* if ( mwt == 0 ) */


  return BPM_SUCCESS;
}
	     
