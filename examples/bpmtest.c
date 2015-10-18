/**
   libbpm-test
   
   Author : Bino Maiheu, University College London (2007)

   Little program to test some routines in libbpm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bpm/bpm_interface.h>
//#include <bpm/bpm_process.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_simulation.h>
#include <bpm/bpm_rf.h>
#include <bpm/bpm_nr.h>

int main( int argc, char **argv ) {

  double xval[10] = {  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0 };
  double yval[10] = { 13.4, 25.6, 32.4, 46.7, 54.5, 62.9, 76.8, 88.2, 95.6, 102.3 };
  double ysig[10] = {  0.3,  0.3,  0.2,  0.5,  0.2,  0.1,  0.7,  0.3,  0.3,   0.1 };
  double a,b,siga,sigb,chi2,q;
  int i;

  /*
  // some rf setup stuff
  rf_setup( 10005, 20.*GHz );
  printf( "rf_nsamples = %d\n", rf_nsamples );
  fflush( stdout );
  */

/*
  for ( i=0; i<1000; i++ ) {
    printf( "rnd[%d] = %f\n", i, nr_rangauss(4.5,1.2) );
    fflush( stdout );
  }
*/



  // test fitting routine...
/*
  nr_fit( xval, yval, 10, NULL, 0, &a, &b, &siga, &sigb, &chi2, &q );
  printf( "fitresults : y = %f + %f x \n", a, b );
  printf( "  siga = %f, sigb = %f \n", siga, sigb );
  printf( "  chi2/ndf = %f\n", chi2 / 8 );
  printf( "  q factor = %f\n", q );
  fflush( stdout );
*/
  // reading some waveforms, calculating some properties :)
/*
  FILE *fp = fopen( "examples/data/run-1456-x10.asc", "r" );
  char buffer[8192], *p;
  int *wf = NULL;
  int evnum, dum, ns;
  double t0, ped, rms;

  double **fft = NULL;

  printf( "%d %d \n", sizeof(double*), sizeof(double) ); 
  fflush( stdout );
  
  while( fgets( buffer, 8192, fp ) ) {
    
    p = strtok( buffer, " " ); evnum = (int) atoi( p );
    p = strtok( NULL, " " ); dum = (int) atoi( p );
    p = strtok( NULL, " " ); dum = (int) atoi( p );
    p = strtok( NULL, " " ); ns = (int) atoi( p );
    
    // alloc wf if for first time...
    if( ! wf ) { wf = (int*) calloc( ns, sizeof(int) ); }
    if( ! fft ) {
      fft = (double**) calloc( ns, sizeof(double*) );
      for(i=0;i<ns;i++) {
	fft[i] = (double*) calloc( 2, sizeof(double) );
      }
    }

    for( i = 0; i < ns ; i++ ) {
      p = strtok( NULL, " " ); wf[i] = (int) atoi( p );
    }
    // get pedestal and t0
    get_pedestal( wf, ns, 20, &ped, &rms );
    if( get_t0( wf, ns, &t0 ) == BPM_SUCCESS ) {
      //printf( "%d %f %f %f \n", evnum, ped, rms, t0 );
      //fflush( stdout );
    }

    // do fft ?
    fft_waveform( wf, ns, fft );
    printf( "%d ", evnum );
    for ( i=0; i<ns; i++ ) printf( "%f ", sqrt( SQR( fft[i][Re] ) + SQR( fft[i][Im] ) ) );
    printf( "\n" );
    fflush( stdout );
    


  }  
  free( wf );
  for( i=0; i<ns; i++ ) free( fft[i] );
  free( fft );
*/
/*
  int wav[256];
  int ev = 0;

  double gen_amp, gen_tdecay, gen_phase, gen_freq;
  double amp, phase, freq, tdecay;

  double **ddc, *w;
  double t0 = 0.2*usecond, t, f;
  double ped, rms, A, C;
  double fwhm = 0.;

  int n1 = 0, n2 = 128;
  int itest;


  ddc = alloc_complex_wave_double( 256 );
  w  = alloc_simple_wave_double( 256 ); 
  nr_seed(12345);

  for (ev=1; ev<=1; ev++ ) {

    gen_amp    = nr_rangauss( 1000., 750. );
    gen_phase  = nr_rangauss( PI/2., PI/2. ); 
    gen_freq   = nr_rangauss( 20.0*MHz, 5.0*MHz );
    gen_tdecay = nr_rangauss( 0.2*usecond, 0.05*usecond );

    //simple_wave( gen_amp, gen_phase, t0, gen_freq, gen_tdecay, 
	//	 pow(2.,13.), 3., 0., 119.*MHz, 14, wav, 256 );


    fft_waveform( wav, 256, ddc );
    fit_fft( ddc, 256, 119.*MHz, &freq, &tdecay, &A, &C );

    for ( i=0; i<128; i++ ) {
      f = (double) i*119.0*MHz/256;
      //sample_to_freq( 119.0*MHz, 256, i, &f );
      
	freq_to_sample( 119.0*MHz, 256, f, &itest );
	printf( "i = %d, itest = %d, freq = %f \n", i, itest, f );
      
      fwhm = 1./PI/tdecay;
      w[i] = A / ( SQR( f - freq ) + SQR( fwhm/2. ) ) + C;
    }

    printf( "%d ", ev );
    for ( i=0; i<128; i++ ) printf( "%f ", SQR(ddc[i][Re])+SQR(ddc[i][Im]) ); printf( "\n" );
    printf( "%d ", ev );
    for ( i=0; i<128; i++ ) printf( "%f ", w[i] ); printf( "\n" );
    fflush( stdout );
    
    fprintf( stderr, "*** gen_freq = %f, gen_tdecay = %f, freq = %f, tdecay = %f, rat = %f \n",
	     gen_freq, gen_tdecay, freq, tdecay, tdecay / gen_tdecay );
    
    
    get_pedestal( wav, 256, 20, &ped, &rms );

    // fit the waveform

    fit_waveform( wav, 256, 0.2*usecond, 119.0*MHz, 20.0*MHz, 0.2*usecond, 1000., 1.5,
		  &freq, &tdecay, &amp, &phase );

    printf( "gen_amp = %f\n", gen_amp );
    printf( "gen_phase = %f\n", gen_phase );
    printf( "fitted pars \n" );
    printf( "freq   = %f\n", freq/MHz );
    printf( "tdecay = %f\n", tdecay/usecond);
    printf( "amp    = %f\n", amp );
    printf( "phase  = %f\n", phase );
    fflush( stdout );
    
    for ( i=0; i<256; i++ ) {
      sample_to_time( 119.0*MHz, 256, i, &t );
      if ( t >= t0 ) {
	w[i] = amp*exp(-(t-t0)/tdecay)*sin( 2.*PI*freq*(t-t0)+phase) + ped;
      } else {
	w[i] = ped;
      }
    }
    printf( "%d ", ev );
    for ( i=0; i<256; i++ ) printf( "%d ", wav[i] ); printf( "\n" );
    printf( "%d ", ev );
    for ( i=0; i<256; i++ ) printf( "%f ", w[i] ); printf( "\n" );
    fflush( stdout );

    printf( "%d %f %f %f %f %f %f\n", ev, gen_amp, gen_phase, amp, phase, freq, tdecay ); 
    fflush( stdout );
  
    downmix_waveform( wav, 256, 119.*MHz, 20.0*MHz, 0.2*usecond, ddc );

    ddc_waveform( wav, 256, 14, 119.0*MHz, 0.2*usecond, 
                  20.0*MHz, 0.2*usecond, 5.0*MHz, 0.001, ddc );
    
    for ( i=0; i<256; i++ ) printf( "%d %d %d %f %f \n", ev, i, wav[i],
				    sqrt( SQR(ddc[i][Re]) + SQR(ddc[i][Im])),
				    atan2( ddc[i][Im], ddc[i][Re] ) ); 
    fflush( stdout );
    

    ddc_sample_waveform( wav, 256, 14, 119.*MHz, 
			 0.2*usecond, 100.*nsecond, // t0, t0Offset
			 20.0*MHz, 0.2*usecond, // freq, tdecay
			 5.0*MHz, 0.001, // filtBW, epsFilt
			 &amp, &phase );
     }

  free_simple_wave_double( w );
  free_complex_wave_double( ddc, 256 );
*/
  return 0;
}
