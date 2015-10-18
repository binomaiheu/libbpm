/*
  libbpm - BPM signal processing/simulation library
  Copyright (C) 2006-07 Bino Maiheu (bino@hep.ucl.ac.uk)
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA 
*/

/** 
    @defgroup dsp Digital Signal Processing Routines
    @anchor dsp
    
    This module contains the definitions for the digital signal
    processing routines for libbpm. 

    @section dsp_filtering The digital filtering routines
    
    @subsection dsp_filtering_usage General usage

    Setup a filter using the create_filter() routine.  
    @code
    filter_t *filter = create_filter( "the_filter", RESONATOR | , 0,
                                       nsamples, 40.*kHz, 8.*kHz, 0., 200. );
    @endcode
    
    The arguments the filter expects is a name for the filter (just for esthetic purposes
    when printing the filter), the filter options, which are explained below, the order
    of the filter, where it is meaning full (e.g. Butterworth, Bessel, Chebyshev). Then
    it needs the number of samples in the waveforms which will be filtered by this filter,
    the sampling frequency and one (optionally two) frequency parameter. For lowpass/highpass
    filters and the resonater, only the first frequency defines respectively the 
    -3dB frequency level for the low/high pass and the resonance frequency for the 
    resonator (the witdh is defined by the Q value in this case). For bandpass/stop filters
    the two frequencies are required and define the -3dB level which defines the bandwidth
    of the filter, with f1 being the lower end frequency and f2 the higher end.

    The implemented filters are :
    - BESSEL      : Bessel IIR filter
    - BUTTERWORTH : Butterwordth IIR filter
    - CHEBYSHEV   : Chebyshev IIR filter
    - RESONATOR   : Resonators
    - GAUSSIAN    : Non-causal Gaussian FIR filter    

    The IIR Bessel, Butterworth and Chebyshev filters can be normalised as lowpass
    (option LOWPASS) which is the default, highpass (option HIGHPASS), bandstop 
    (option BANDSTOP) or bandpass (option BANDPASS) filters. They are designed with 
    poles and zeros in the s plane that are transformed to the z plane either by bilinear 
    z transform (option BILINEAR_Z_TRANSFORM) or matched z transform (option MATCHED_Z_TRANSFORM).
    Just "OR" the options together to setup the filter, e.g. :

    @code
    filter_t *filter = create_filter( "lp", BESSEL | HIGHPASS | MATCHED_Z_TRANSFORM, 0,
                                       ns, 40.*kHz, 8.*kHz, 0., 200. );
    @endcode

    The resonators are designed directly with their 2 poles and 2 zeros in the z plane and
    can be normalised either as BANDPASS (default), BANDSTOP (or NOTCH) or ALLPASS
    resonators.

    The last argument to the create_filter() routine is a parameter which can optionally be
    given to the filter. It depends on the filter chosen, currently the parameter has
    meaning for the following filters :

    - BESSEL : the parameter defines the ripple in dB, has to be negative !
    - RESONATOR : the parameter gives the Q value of the resonator, if you want to have
                  a pure oscillator (so infinite Q), then set the parameter to a negative
		  number or zero.
    - GAUSSIAN : the filter cut-off parameter, or the fraction of the gaussian convolution 
                 function below which it is set to 0. (default is 0.001)
		  
    The filter coefficients for the difference equation are calculated and checked
    for consistency, upon which they are stored in the filter structure. Once this
    is done and the filter is setup, application to various waveforms is fairly 
    straightforward. Note that you only have to define your filter once during 
    initialisation. Once setup, it can be used to filter any number of waveforms of
    the same type.

    @code
    apply_filter( filter, wave );
    @endcode

    
    To get an impulse response from the filter into the secified waveform, where
    the impulse is given at sample 1000, the following routine is implemented. 

    @code
    filter_impulse_response( filter, wave, 1000 );
    @endcode

    This routine creates an impulse function (zero everywhere, except at the sample 
    you enter, where it's value is 1) and puts it through the filter. The FFT of this
    impulse response gives you the filter characteristic in frequency domain. Also you
    can check the filter's response to a step function, it's so-called step response :

    @code
    filter_step_response( filter, wave, 1000 )
    @endcode

    The step response is defined as the response of the filter to an input function
    which is zero at the beginning and 1 for samples >= the sample you specify.


    @subsection dsp_filtering_bbc The Bessel, Butterworth and Chebyshev filters

    
    @subsection dsp_filtering_resonator The Resonator filter


    @subsection dsp_filtering_gaus The gaussian filter
    The gaussian filter is implemented as a FIR convolution with both causal and anti-causal
    coefficients. Note that the frequency given is treated as the -3dB level for the 
    gaussian. There is an option to restore the definition for bandwith which was 
    used in early ESA processing, being the gaussian sigma, use GAUSSIAN_SIGMA_BW.


    @section dsp_ddc The Digital Downconversion Algorithm (DDC)

    The digital downconversion routine was developped to process digitised BPM waveforms
    and to retreive their position and amplitude. It basically implements an RF mixer in 
    software. You need to supply it with the doublewf_t holding the waveform to mix down
    and the frequency for the software LO. Also you need to give a pointer to a
    low-pass filter in order to filter out the resulting double frequency component from
    the downmixing. The routine

    @code
    int ddc( doublewf_t *w, double f, filter_t *filter, complexwf_t *dcw );
    @endcode

    returns then the complex DC waveform (dcw), where it's amplitude and phase can then
    be used in further calculations for beam position and slope in the BPM. We recommend
    the usage of a GAUSSIAN low-pass filter for the double frequency filtering as this 
    shows the best phase behaviour combined with linearity (see create_filter()).

    For fast execution, the DDC routine comes with a buffer which it only allocates once
    by doing 

    @code
    ddc_initialise();
    @endcode
    
    This buffer is used in the filtering routine, you can clean up after the execution
    of the buffer by having

    @code
    ddc_cleanup();
    @endcode

    @section dsp_fft Discrete (Fast) Fourier Transforms

    The FFT routines in the dsp section of libbpm are based upon the General Purpose FFT
    Package by Takuya OOURA, 1996-2001, see http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
    More specifically on it's  split-radix fast version (fftsg). These set of routines
    needs a buffer for bitswapping an a buffer to store a table with sin and cos values
    so they needn't be calculated for every FFT. The routine

    @code
    fft_initialise( int ns )
    @endcode

    intialises the buffers for waveforms of a certain sample length ns. Note that ns has to
    be a power of 2. You can clear the FFT buffers by issuing

    @code
    fft_cleanup( );
    @endcode

    Then two wrapper routines are implemented which take doublewf_t and complexwf_t data.

    @subsection dsp_fft_complex Complex Discrete Fourier Transform
    The first one is

    @code
    int complexfft( complexwf_t *z, int fft_mode );
    @endcode

    which takes a complex waveform and performs an FFT in place. The fft_mode argument
    can be either 
    - FFT_FORWARD  : forward discrete Fourier transform (plus-sign)

    @f[ X[k] = \sum_{j=0}^{n-1} x[j]*\exp(2*\pi*i*j*k/n), 0<=k<n @f]

    - FFT_BACKWARD : backward discrete Fourier transform (minus-sign)

    @f[ X[k] = \sum_{j=0}^{n-1} x[j]*\exp(-2*\pi*i*j*k/n), 0<=k<n @f]
    
    Note the backward and forward FFT's have a factor of n inbetween them, so to get the
    orginal wf back after applying both the backward and the forward FFT, you need
    to divdide by the number of samples z->n. 
    
    @subsection dsp_fft_real Real Discrete Fourier Transform
    The second routine implements the real discrete Fourier transform when
    having FFT_FORWARD and the other way around when having FFT_BACKWARD.

    @code
    int realfft( doublewf_t *y, int fft_mode, complexwf_t *z );
    @endcode

    So for FFT_FORWARD

    @f[ Re( X[k] ) = \sum_{j=0}^{n-1} a[j]*\cos(2*\pi*j*k/n), 0<=k<=n/2 @f]
    @f[ Im( X[k] ) = \sum_{j=0}^{n-1} a[j]*\sin(2*\pi*j*k/n), 0<k<n/2   @f]

    and FFT_BACKWARD takes the input frmo the first half (n/2) of the complexwf_t
    and FFTs it, expanding to a doublewf_t of length n. 
    
    @f[ X[k] = \frac{(Re(x[0]) + Re(x[n/2])*\cos(\pi*k) )}{2} + 
               \sum_{j=1}^{n/2-1} Re(x[j])*\cos(2*\pi*j*k/n) + 
               \sum_{j=1}^{n/2-1} Im(x[j])*\sin(2*\pi*j*k/n), 0<=k<n @f]
    
    @subsection dsp_fft_ref Reference for FFT routines
    - Masatake MORI, Makoto NATORI, Tatuo TORII: Suchikeisan, 
      Iwanamikouzajyouhoukagaku18, Iwanami, 1982 (Japanese)
    - Henri J. Nussbaumer: Fast Fourier Transform and Convolution 
      Algorithms, Springer Verlag, 1982
    - C. S. Burrus, Notes on the FFT (with large FFT paper list)
      http://www-dsp.rice.edu/research/fft/fftnote.asc
      
    @subsection dsp_fft_copy Copyright statement for FFT routines
    Copyright(C) 1996-2001 Takuya OOURA
    email: ooura@mmm.t.u-tokyo.ac.jp
    download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
    You may use, copy, modify this code for any purpose and 
    without fee. You may distribute this ORIGINAL package.


    @section dsp_examples DSP example program
    There is an example program, which can be found in the examples directory under dsp.
    It shows how to work with the filtering and the DDC routines...

    @code
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    #include <math.h>

    #include <iostream>
    
    #include <TROOT.h>
    #include <TFile.h>
    #include <TTree.h>
    
    #include <bpm/bpm_process.h>
    #include <bpm/bpm_units.h>
    #include <bpm/bpm_simulation.h>
    #include <bpm/bpm_nr.h>
    #include <bpm/bpm_rf.h>
    #include <bpm/bpm_dsp.h>
    #include <bpm/bpm_wf.h>
    
    using namespace std;

    int main( int argc, char **argv ) {
    
      cout << "Welcome to the libbpm DSP sandbox" << endl;
    
      int ns    = 256;
      double fs = 119.*MHz;
      
      doublewf_t *w = doublewf( ns, fs );
      doublewf_t *s = doublewf_sample_series( ns, fs );
      
      doublewf_t *ddc_amp   = doublewf( ns, fs );
      doublewf_t *ddc_phase = doublewf( ns, fs );
      
      // setup the root trees...
      TFile *rootfile = new TFile( "dsp.root", "recreate" );
      TTree *roottree = new TTree( "dsp", "libbpm dsp tests" );
      
      int evt;
      double amp, phase;
      double gen_amp, gen_phase;
      
      // setup the branches in the tree
      roottree->Branch( "evt",       &evt,          "evt/I"            );
      roottree->Branch( "wf",        w->wf,         "wf[256]/D"        );
      roottree->Branch( "s",         s->wf,         "s[256]/D"         );
      roottree->Branch( "gen_amp",   &gen_amp,      "gen_amp/D"        );
      roottree->Branch( "gen_phase", &gen_phase,    "gen_phase/D"      );
      roottree->Branch( "ddc_amp",   ddc_amp->wf,   "ddc_amp[256]/D"   );
      roottree->Branch( "ddc_phase", ddc_phase->wf, "ddc_phase[256]/D" );
      
      complexwf_t *ddcwf = complexwf( ns, fs ); 

      filter_t *gauss  = create_filter( "gauss", GAUSSIAN,0,ns,fs,6.*MHz,0.,0.001);
      filter_t *butter = create_filter( "butter", BUTTERWORTH | LOWPASS,4,ns,fs,6.*MHz,0.,0.);
      filter_t *bessel = create_filter( "bessel", BESSEL | LOWPASS,4,ns, fs,6.*MHz,0., 0.);
      filter_t *cheby  = create_filter( "cheby",  CHEBYSHEV | LOWPASS,4,ns,fs,6.*MHz,0.,-10.);
      
      // init the DDC
      ddc_initialise( ns, fs );
      
      for ( evt = 1; evt<=1000; evt++ ) {
      
        // Make the waveform
	gen_amp   = (double) evt * 10.;
	gen_phase = PI / (double) evt;
	
	// reset the w to 0... quite important :D
	doublewf_reset( w );
	
	doublewf_add_dcywave( w, gen_amp, gen_phase, 21.4*MHz, 0.15*usec, 0.2*usec, 0. );
	
	// do the DDC :)
	if ( ddc( w, 21.4*MHz, gauss, ddcwf ) ) return 1;
	
	// want to try differen filters ?
	//if ( ddc( w, 21.4*MHz, butter, ddcwf ) ) return 1;
	//if ( ddc( w, 21.4*MHz, bessel, ddcwf ) ) return 1;
	//if ( ddc( w, 21.4*MHz, cheby, ddcwf ) ) return 1;
	
	
	// get amplitude and phase from complex wf
	complexwf_getamp( ddc_amp, ddcwf );
	complexwf_getphase( ddc_phase, ddcwf );
	
	// fill the tree...
	roottree->Fill();
	
	if ( evt % 100 == 0 ) cout << "Simulated " << evt << " events." << endl;
      }    

      // clear the DDC memory buffers
      ddc_cleanup();
      
      rootfile->Write();
      rootfile->Close();
      
      delete_filter( gauss );
      delete_filter( butter );
      delete_filter( bessel);
      delete_filter( cheby );
      
      complexwf_delete( ddcwf );
      
      doublewf_delete( w );
      doublewf_delete( s );
      doublewf_delete( ddc_amp );
      doublewf_delete( ddc_phase );
      
      return 0;
    }
    @endcode
*/

/**
    @file     
    @ingroup dsp
    @brief   libbpm digital signal processing routines
*/

/** @addtogroup dsp */
/** @{ */

#ifndef BPMDSP_H__
#define BPMDSP_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bpm/bpm_defs.h"
#include "bpm/bpm_messages.h"
#include "bpm/bpm_nr.h"
#include "bpm/bpm_wf.h"

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */
#define BESSEL                 0x0000001 /**< Bitmask for Bessel filter */
#define BUTTERWORTH            0x0000002 /**< Bitmask for Butterworth filter */
#define CHEBYSHEV              0x0000004 /**< Bitmask for Chebyshev filter */
#define RAISEDCOSINE           0x0000008 /**< Bitmask for Raised Cosine filter */
#define RESONATOR              0x0000010 /**< Bitmask for Resonator filter */
#define GAUSSIAN               0x0000020 /**< Bitmask for Gaussian filter */

#define BILINEAR_Z_TRANSFORM   0x0000100 /**< Get z poles via bilinear z transform from s plane */
#define MATCHED_Z_TRANSFORM    0x0000200 /**< Get z poles via matches z transform from s plane */
#define NO_PREWARP             0x0000400 /**< Don't do the prewarp correction */
#define CAUSAL                 0x0000800 /**< Filter is purely causal (only depends on past ) */
#define ANTICAUSAL             0x0001000 /**< .... purely anticausal (only depends on future) */
#define NONCAUSAL              0x0001800 /**< Filter is both causal and acausal */
#define GAUSSIAN_SIGMA_BW      0x0002000 /**< Gaussian sigma bandwidth in stead of -3 dB (def)*/

#define LOWPASS                0x0010000 /**< Normalise filter as lowpass */
#define HIGHPASS               0x0020000 /**< Normalise filter as highpass */
#define BANDPASS               0x0040000 /**< Normalise filter as bandpass */
#define BANDSTOP               0x0080000 /**< Normalise filter as bandstop */
#define NOTCH                  0x0080000 /**< Normalise filter as notch filter (=bandstop) */
#define ALLPASS                0x0100000 /**< Normalise filter as allpass ( resonator ) */

#define FIR                    0x0200000 /**< Filter is of FIR type */
#define IIR                    0x0400000 /**< Filter is of IIR type */

#define MAXORDER                   10  /**< Maximum filter order */
#define MAXPZ                      50  /**< Maximum number of poles and zeros >2*MAXORDER */
#define FILT_EPS              1.0e-10  /**< A small number used in bpmdsp */
#define MAX_RESONATOR_ITER         50  /**< Maximum iterations in resonator poles calculation */

#define FFT_FORWARD                 0  /**< Perform FFT from time -> frequency */
#define FFT_BACKWARD                1  /**< Perform FFT from frequency -> time */

/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

/**
 * The filter representation in the complex plane (poles/zeros).
 */
typedef struct {
  int        npoles;       /**< The number of filter poles */
  int        nzeros;       /**< The number of filter zeros */
  complex_t  pole[MAXPZ];  /**< Array of the filter's complex poles */
  complex_t  zero[MAXPZ];  /**< Array of the filter's complex zeros */
} filterrep_t;

/**
 * The filter structure.
 */
typedef struct {
  char          name[80];      /**< The filter's name */

  unsigned int  options;       /**< type and option bits for filter */
  int           order;         /**< filter order                    */
  
  double        fs;            /**< sampling frequency */
  double        f1;            /**< first frequency ( left edge for bandpass/stop ) */
  double        f2;            /**< right edge for bandpass/stop ( undef for low/highpass ) */

  double        alpha1;        /**< rescaled f1 */
  double        alpha2;        /**< rescaled f2 */

  double        w_alpha1;      /**< warped alpha1 */
  double        w_alpha2;      /**< warped alpha2 */
  
  double        cheb_ripple;   /**< ripple for chebyshev filters */
  double        Q;             /**< Q factor for resonators */
  double        gauss_cutoff;  /**< gaussian filter cutoff parameter */

  complex_t     dc_gain;       /**< Complex DC gain of the filter               */
  complex_t     fc_gain;       /**< Complex Center frequency gain of filter     */
  complex_t     hf_gain;       /**< Complex High frequency (fNy) gain of filter */
  double        gain;          /**< Actual Filter gain */
  
  filterrep_t  *cplane;        /**< pointer to complex filter representation, poles and zeros */

  int           nxc;            /**< number of x coefficients */
  double        xc[MAXPZ+1];    /**< pointer to array of x coefficients */

  int           nxc_ac;         /**< number of anti-causal x coefficients */
  double        xc_ac[MAXPZ+1]; /**< pointer to array of anti-causal x coefficients */

  int           nyc;            /**< number of y coefficients (for IIR filters) */
  double        yc[MAXPZ+1];    /**< pointer to array of y coefficients */

  int           nyc_ac;         /**< number of anti-causal y coefficients (for IIR filters) */
  double        yc_ac[MAXPZ+1]; /**< pointer to array of anti-causal y coefficients */

  double        xv[MAXPZ+1];     /**< filter x buffer, used in apply_filter */
  double        xv_ac[MAXPZ+1];  /**< filter x buffer, used in apply_filter */

  double        yv[MAXPZ+1];     /**< filter y buffer, used in apply_filter */
  double        yv_ac[MAXPZ+1];  /**< filter y buffer, used in apply_filter */

  int           ns;            /**< number of samples of waveforms to be filtered */
  double       *wfbuffer;      /**< waveform buffer for filter computations, allocated once ! */
  
} filter_t;

  
/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */

  /**
     Creates the filter.
     @param name    a name for the filter
     @param options filter specification and options bitword 
     @param order   filter order
     @param ns      number of samples of the waveforms
     @param fs      sampling frequency
     @param f1      first frequency
     @param f2      optional second frequency ( bandpass/bandstop )
     @param par     optional parameter
                    - for chebyshev : ripple in dB
                    - for resonator : Q factor
     @return A pointer to the created filter structure, memory is allocated 
             on the heap inside this routine, the user has to take of deleting
             it using delete_filter().
  */
  EXTERN filter_t* create_filter( char name[], unsigned int options, int order,
				  int ns, double fs, double f1, double f2, double par );

  /**
     Apply the filter to the given waveform.
     Note that the filter is applied in place, the user has to make a copy
     of the waveform if he/she wants to keep the original before applying the
     filter. The number of samples in the waveform has to be set in advance
     when creating the filter, it is stored in the filter structure (f->ns).

     @param  f   pointer to a filter that was created using create_filter
     @param  wf  an array containing the waveform to be filtered
     @return BPM_SUCCESS upon success and BPM_FAILURE upon failure
  */
  EXTERN int apply_filter( filter_t *f, doublewf_t *w );


  /**
     Prints the filter to the given file pointer.
     @param of the filepointer, use "stdout" to print to the terminal
     @param f  the filter to be printed
     @return void
  */  
  EXTERN void print_filter( FILE *of, filter_t *f );

  /**
     Clears the memory that was allocated on the heap for the filter f.
     @param f a pointer to the filter
     @return void
   */
  EXTERN void delete_filter( filter_t *f );


  /**
     This routine fills the given wf with the step response of the filter. The step 
     response is defined as wf[i] = 0. for i < itrig and wf[i] = 1. for i >= itrig.

     @param f     a pointer to the filter to use
     @param wf    pointer to a waveform which will be overwritten with the step response
     @param itrig the sample number in the waveform which will have the step
     @return BPM_SUCCESS upon succes and BPM_FAILURE upon failure
  */
  EXTERN int filter_step_response( filter_t *f, doublewf_t *w, int itrig );

  /**
     This routine fills the given wf with the impulse response of the filter. The 
     impulse response is defined as wf[i] = 1. for i == itrig and wf[i] = 0. elsewhere.

     @param f     a pointer to the filter to use
     @param wf    pointer to a waveform which will be overwritten with the impulse response
     @param itrig the sample number in the waveform which will have the impulse
     @return BPM_SUCCESS upon succes and BPM_FAILURE upon failure
  */
  EXTERN int filter_impulse_response( filter_t *f, doublewf_t *w, int itrig );
  


  /**
     This routine returns a pointer to a filter representation filterrep_t in the 
     s plane for Butterworth, Chebyshev and Bessel filters. It need an initialised
     filter structure which has the filter type and the order set. Memory
     is allocated for this routine on the heap, so the user is responsible to 
     delete this memory using free(). 

     @param f the initialised filter with the correct options in f->options
     @return the filter representation in the s plane
  */
  EXTERN filterrep_t* create_splane_representation( filter_t *f );


  /**
     This routine returns a pointer to a filter representation filterrep_t in the 
     z plane for resonance filters. It needs an initialised filter structure which has 
     the filter type and the Q factor set. Memory is allocated for this routine on the
     heap, so the user is responsible to delete this memory using free(). 

     @param f the initialised filter with the correct options in f->options
     @return the filter representation in the z plane
  */
  EXTERN filterrep_t* create_resonator_representation( filter_t *f );

  /**
     This routine transforms the poles and zeros for Bessel, Chebyshev and Butterworth
     filters to the z plane either via matched z transform or bilinear z transform.
     This is set in f->options. Memory is allocated for this routine on the
     heap, so the user is responsible to delete this memory using free().

     @param f the filter, needs the options from it to check how to transform
     @param s filter s plane poles and zeros
     @return a pointer to the z plane representation
  */
  EXTERN filterrep_t* zplane_transform( filter_t *f, filterrep_t *s );


  /**
     Prints the filter representation in terms of poles and zeros to the filepointer.
     @param of the filepointer, use "stdout" to print to the terminal
     @param r  the filter representation to be printed
     @return void
   */
  EXTERN void print_filter_representation( FILE* of, filterrep_t *r );

  
  /**
     Normalises the Butterworth, Chebyshev or Bessel filters to be Bandpass/stop or Low/Highpass
     @param f  the filter
     @param s  the filter's representation in the s plane
     @return BPM_SUCCESS upon success or BPM_FAILURE upon failure.
   */
  EXTERN int normalise_filter( filter_t *f, filterrep_t *s );


  /**
     Calculates the filter coefficients from the z plane representation for
     Butterworth, Chebyshev, Bessel and Resonators. Before this routine is called, 
     one has to make sure that the member cplane, which holds a pointer to the 
     filter's representation in the complex plane is set. This routine than calculates 
     the filter coefficients and stores them in 
     f->xc ( coefficients of x[n], x[n-1], x[n-2]...) and f->yc ( coefficients of 
     y[n-1], y[n-2], y[n-3], ... in case of IIR filters ). 

     @param f the filter, having it's f->cplane member set to the z plan representation
     @return BPM_SUCCESS upon success or BPM_FAILURE upon failure.
  */
  EXTERN int calculate_filter_coefficients( filter_t *f );
  

  /** 
      Calculates the gaussian filter coefficients from the original gaussian filter
      implementation in the digital downconversion algortithm in Yury's code.
      Note that this filter is implemented as a FIR non-causal filter. 
      @param f the filter structure with the coefficients to fill
      @return BPM_SUCCESS upon success or BPM_FAILURE upon failure.
  */
  EXTERN int gaussian_filter_coeffs( filter_t *f );


  /**
     Helper routine to expand a complex polynomial from a set of zeros.
     @param w array of complex zeros for the polynomial
     @param n nunber of zeros
     @param a array of coeffiecients for the polynomial that is returned
     @return BPM_SUCCESS upon success or BPM_FAILURE upon failure.
  */
  EXTERN int _expand_complex_polynomial( complex_t *w, int n, complex_t *a );

  /**
     Helper routine to evaluate a complex polynomial for value z
     @param a array of coeffiecients for the polynomial that is returned
     @param n number of zeros
     @param z the value for which to evalute the polynomial
     @return the value of the polynomial for z ( complex_t )
  */
  EXTERN complex_t _eval_complex_polynomial( complex_t *a, int n, complex_t z );


  // -- DDC routines -----------------------------------------------------

  /**
     Initialises and allocates memory for the DDC buffers with the correct 
     number of samples and sampling frequency 
     @param ns Nuber of samples in waveforms to be processed 
     @param fs The sampling frequency of the waveforms
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int ddc_initialise( int ns, double fs );

  /**
     Clears up and frees the buffer memory for the ddc routines */
  EXTERN void ddc_cleanup( void );
  
  /**
     Do a digital downconversion on the waveform f. The routine
     returns a complex DC waveform "wdc". If the buffer arguments
     are NULL pointers, the DDC routine will use an internal buffer. 
     This is a good option when all the BPMs in the system have
     the same sampling frequency and number of samples.
     @param w The waveform of doubles to process
     @param f The frequency of the digital local oscillator
     @param filter The lowpass filter to get rid of the 2omega component
     @param dcw The complex DC waveform
     @param bufre The real ddc buffer
     @param bufim The imaginary ddc buffer
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
   */
  int ddc( doublewf_t *w, double f, filter_t *filter, complexwf_t *dcw,
	   doublewf_t *bufre, doublewf_t *bufim );


  // -- DFT (FFT) routines ---------------------------------------------------

  /**
     Regenerates the sin/cos tables that are needed for the fast DFT algorithm. */
  EXTERN int  fft_gen_tables( void );

  /**
     This one initialised the FFT buffers, checks whether they are large enough for the
     given number of samples and frees and re-allocates memory where necessary 
     @param ns The number of samples in the waveforms to be transformed
     @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int  fft_initialise( int ns );
  
  /**
     This routine frees up the memory used by the FFT buffers
  */
  EXTERN void fft_cleanup( void );
  
  /** 
      Executes a complex fast fourier transform in line. See the reference
      guide for details.
      @param z The complex waveform to transform (original waveform is destroyed)
      Note that the number of samples need to be a power of 2.
      @param fft_mode Specifies whether to do the forward or backward transform
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int  complexfft( complexwf_t *z, int fft_mode ); 

  /**
     Executes a real fast fourier transform, between the real waveform y and the 
     complex waveform z. See documentation for further explanation.
     @param y Pointer to the real wavefrom
     @param fft_mode Specifies whether to do the forward or backward transform
     @param z Pointer to the complex waveform
     @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int  realfft( doublewf_t *y, int fft_mode, complexwf_t *z );

  /**
     Normalises the phase, to the interval [0,2pi[ 
     @param phase Pointer to the phase value to normalise
  */
  EXTERN void norm_phase( double *phase );

#ifdef __cplusplus
}
#endif

#endif /* #ifndef __BPM_FILTER_H_LOADED */
/*@}*/
/* ================================ end of file ============================= */

