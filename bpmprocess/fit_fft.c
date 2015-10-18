/**
   @file
   @ingroup processing
*/
#include <stdio.h>
#include <stdlib.h>
#include <bpm/bpm_nr.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

#define FIT_MAX_ITER     5000  // max iteration
#define FIT_WINDOW_FACTOR   3  // number of sigma to fit withing (integer)

/*
  p[0] = amplitude
  p[1] = freq
  p[2] = tdecay
  p[3] = offset
*/

void fcnlorjac( double *p, double *ljac, int np, int ns, void *a ) {

  int i,j;
  double denom;
  double f;

  for ( i=j=0; i<ns; i++ ) {
    // convert sample nr to freq
    f = (((double*)a)[0] + (double)i + 0.5 ) / ((double*)a)[1] * ((double*)a)[2];
    denom = SQR( f - p[1] ) + SQR( p[2] / 2. );

    if ( denom == 0. ) {
      ljac[j++] = 0.;
      ljac[j++] = 0.;
      ljac[j++] = 0.;
      ljac[j++] = 1.;
    } else {
      ljac[j++] =  1. / denom;
      ljac[j++] =  2. * p[0] * (f-p[1]) / SQR( denom );
      ljac[j++] = -0.5 * p[0] * p[2] / SQR( denom );
      ljac[j++] =  1.;
    }

  }

  return;
}

void fcnlor( double *p, double *lor, int np, int ns, void *a ) {

  int i;
  double f;
  
  // We'll do this fit in terms of samples and inverse samples, convert back
  // in fit_fft to time and frequency in usec and MHz
  // a[0] = (double) n1
  // a[1] = ns;
  // a[2] = fs;

  for ( i=0; i<ns; i++ ) {
    // convert sample nr to freq
    f = (((double*)a)[0] + (double) i + 0.5 ) / ((double*)a)[1] * ((double*)a)[2];
    lor[i] = p[0] / ( SQR( f - p[1] ) + SQR( p[2] / 2. ) ) + p[3];
  }

  return;
}

// ----------------------------------------------------------------------------------

int fit_fft_prepare( complexwf_t *ft, int *n1, int *n2, 
		     double *amp, double *freq, double *fwhm ) {

  int i, imax=0;
  double  pw, tmp;

  if ( ! ft || ! amp ) {
    bpm_error( "Invalid pointers in fit_fft_prepare(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // some vary basic initial values
  *amp    = 0; // need this to be 0 to find max !
  *freq   = 20.0 * _MHz__;
  *fwhm   = 10.0 * _MHz__; 
  *n1     =  20;
  *n2     = 100;

  // get maximum peak position in spectrum, only up to nyquist freq (ns/2)
  for ( i=0; i<ft->ns/2; i++ ) {
    pw = c_abs2( ft->wf[i] ); // SQR( fft[i][Re] ) + SQR( fft[i][Im]);
    if ( pw > *amp ) {
      *amp = pw; imax = i;
    }
  }
  *freq = (double) imax / (double) ft->ns * ft->fs;

  // run left and right to get the FWHM initial value
  for ( i=imax; i>0; i-- ) {
    pw = c_abs2( ft->wf[i] ); // SQR( fft[i][Re] ) + SQR( fft[i][Im]);
    if ( pw <= *amp/2. ) break;
  }
  *n1 = i;

  for ( i=imax; i<ft->ns/2 ; i++ ) {
    pw = c_abs2( ft->wf[i] ); //SQR( fft[i][Re] ) + SQR( fft[i][Im]);
    if ( pw <= *amp/2. ) break;
  }
  *n2 = i;

  // set the FWHM initial decay cte
  *fwhm = ( (double)( *n2 - *n1 ) / (double) ft->ns * ft->fs );

  // now let's take 2 x FWHM as fit range...
  *n1 = imax - FIT_WINDOW_FACTOR * ABS(imax - *n1);
  *n2 = imax + FIT_WINDOW_FACTOR * ABS(imax - *n2);

  // safety check...
  if ( *n1 < 0 ) *n1 = 0;
  if ( *n2 > ft->ns/2) *n2=ft->ns/2;

  if ( *n2 <= *n1 ) {
    bpm_error( "Error in fit range ( n2 <= n1 ) in fit_fft_prepare(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  if ( *n2 - *n1 < 5 ) {
    bpm_error( "Error, too few number of samples in fit_fft_prepare(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }


#if 0
  printf( "fftfit prepare : n1=%d, n2=%d, imax=%d, amp=%f, freq=%f, fwhm=%f\n",
	  *n1, *n2, imax, *amp, *freq, *fwhm );
  fflush( stdout );
#endif

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------------

int fit_fft( complexwf_t *ft, double *freq, double *tdecay, double *A, double *C ) {

  const int npar = 4; // number of parameters in fit
  double *pspec, par[npar];
  double opt[4], info[LM_INFO_SZ];
  double a[3], i_amp, i_freq, i_fwhm;
  int i, nfit, n1, n2;

#if 0
  double *err;
#endif

  opt[0]=LM_INIT_MU; opt[1]=1E-15; opt[2]=1E-15; opt[3]=1E-20;

  // init
  *freq   = 0.;
  *tdecay = 0.;

  if ( ! ft ) {
    bpm_error( "Invalid pointer in fit_fft(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  if ( fit_fft_prepare( ft, &n1, &n2, &i_amp, &i_freq, &i_fwhm ) == BPM_FAILURE ) {
    return BPM_FAILURE;
  }

  // prepare the fit stuff
  nfit = n2 - n1 + 1;
  pspec = (double*) calloc( nfit, sizeof(double) ); 
  for ( i=0; i<nfit; i++ ) {
    pspec[i] = c_abs2( ft->wf[i] ); // SQR( fft[n1+i][Re] ) + SQR( fft[n1+i][Im] );
  }

  // intial pars
  par[0] = i_amp;
  par[1] = i_freq;
  par[2] = i_fwhm;
  par[3] = 0.;

  // extra data
  a[0] = (double) n1;
  a[1] = (double) ft->ns;
  a[2] = ft->fs;

  // test jacobian
#if 0
  err = (double*) calloc( nfit, sizeof(double) );
  nr_lmchkjac( fcnlor, fcnlorjac, par, npar, nfit, a, err );
  for (i=0;i<nfit;i++) printf( "err[%d] = %f\n", i, err[i] ); fflush(stdout);
  free(err);
#endif

  // here we go
  nr_lmder( fcnlor, fcnlorjac, par, pspec, npar, nfit, FIT_MAX_ITER, opt, info,
		    NULL, NULL, a );

  *freq = par[1];
  // we fitted the FWHM, this is the inverse of the decay time and lacks a 
  // factor of pi : http://mathworld.wolfram.com/FourierTransformLorentzianFunction.html
  if ( par[2] != 0. ) *tdecay = 1. / par[2] / PI;

  // return the other parameters as well if needed
  if ( A ) *A = par[0];
  if ( C ) *C = par[3];
  
  free( pspec );

  return BPM_SUCCESS;
}

#undef FIT_MAX_ITER
