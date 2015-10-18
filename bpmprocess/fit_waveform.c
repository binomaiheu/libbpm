/**
   @file
   @ingroup processing
*/
#include <bpm/bpm_nr.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_process.h>

#define FIT_MAX_ITER 10000

#define FIT_AMP          0
#define FIT_PHASE        1
#define FIT_FREQ         2
#define FIT_TDECAY       3

#define FIT_T0           0
#define FIT_FS           1

void fcnwfjac( double* par, double *jac, int npars, int ns, void *a ) {
  int i, j;
  double t, sinarg, exparg, cosarg;

  // we trust we have npars = 4 here !!

  for( i=j=0; i<ns; i++ ) {
    sample_to_time( ( (double*)a )[FIT_FS], ns, i, &t );
    if ( t >=  ((double*)a)[FIT_T0] ) {
      sinarg = sin( 2.*PI*par[FIT_FREQ]* ( t - ((double*)a)[FIT_T0] ) + par[FIT_PHASE] );
      cosarg = cos( 2.*PI*par[FIT_FREQ]* ( t - ((double*)a)[FIT_T0] ) + par[FIT_PHASE] );
      exparg = exp( -( t - ((double*)a)[FIT_T0] )  / par[FIT_TDECAY] );
      
      // these have to be in order of appearance :)

      // \partial x / \partial p0 (amp)
      jac[j++] = exparg * sinarg;

      // \partial x / \partial p1 (phase)
      jac[j++] = par[FIT_AMP] * exparg * cosarg; 

      // \partial x / \partial p2 (freq)
      jac[j++] = par[FIT_AMP] * exparg * cosarg * 2.*PI*( t - ((double*)a)[FIT_T0] );

      // \partial x / \partial p3 (tdecay)
      jac[j++] = par[FIT_AMP] * exparg * (t - ((double*)a)[FIT_T0])/SQR( par[FIT_TDECAY]) * sinarg;
      
    } else {
      jac[j++] = 0.;
      jac[j++] = 0.;
      jac[j++] = 0.;
      jac[j++] = 0.;
    }
    
  }
}

/*
   The fitfunction, being simply the waveform, setup for the additional data array 
   xval[0] = t0
   xval[1] = the sampling frequency
*/
void fcnwf( double* par, double *sinwf, int npars, int ns, void *a ) {

  double t;
  int i;

  for( i=0; i<ns; i++ ) {
    sample_to_time( ( (double*)a )[FIT_FS], ns, i, &t );
    if ( t >=  ( (double*)a )[FIT_T0] ) {
      sinwf[i] = par[FIT_AMP] 
	* exp( - ( t-( (double*)a )[FIT_T0]) / par[FIT_TDECAY] ) 
	* sin( 2.*PI*par[FIT_FREQ]*(t-( (double*)a )[FIT_T0]) + par[FIT_PHASE] );
    } else {
      sinwf[i] = 0.;
    }
  }

  return;
}

int fit_waveform( doublewf_t *w, double t0, 
		  double i_freq, double i_tdecay, double i_amp, double i_phase,
		  double *freq, double *tdecay, double *amp, double *phase ) {

  doublewf_t *yval;
  double a[2];
  const int npars = 4;
  double par[npars]; // 4 parameters in fit : amp, phase, freq, tdecay
  double opt[4], info[LM_INFO_SZ];
  int i;

  double err[256];

  opt[0]=LM_INIT_MU; opt[1]=1E-15; opt[2]=1E-15; opt[3]=1E-20;

  if ( ! w ) {
    bpm_error( "Invalid waveform pointer in fit_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // need a double array
  yval = doublewf( w->ns, w->fs );
  if ( ! yval ) {
    bpm_error( "Unable to allocate memory for waveform in fit_waveform(...)",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // preparing
  a[FIT_T0] = t0;
  a[FIT_FS] = w->fs;

  par[FIT_AMP]    = i_amp;
  par[FIT_PHASE]  = i_phase;
  par[FIT_FREQ]   = i_freq;
  par[FIT_TDECAY] = i_tdecay;
  
  // we've calculated the jacobian analytically, so let's use it !
  nr_lmder( fcnwf, fcnwfjac, par, yval->wf, npars, w->ns, FIT_MAX_ITER, opt, info, NULL, NULL, a );

  /*
    TODO :
    - Probably need to put in some box constraints for better behaviour at a certain point, 
      can use the routine nr_lmder_bc for this. 
    - Maybe program in something to replace the matrix inversion based upon the LU algorithm
      in nr_levmar.c. The original version had some LAPACK functionality, but as we
      wanted this standalone, i've cut it out. 
    - Best way to use it probably to take the parameters from previous fit of last pulse and have
      them as initials for the current fit.
  */


  *amp    = par[FIT_AMP];
  *phase  = par[FIT_PHASE];
  *freq   = par[FIT_FREQ];
  *tdecay = par[FIT_TDECAY];

  doublewf_delete( yval );

  return BPM_SUCCESS;
}

#undef FIT_MAX_ITER
#undef FIT_AMP      
#undef FIT_PHASE    
#undef FIT_FREQ     
#undef FIT_TDECAY   
#undef FIT_T0       
#undef FIT_FS

