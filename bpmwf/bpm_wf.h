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
   @defgroup wave Waveform handling routines
   @anchor wave

   This module contains the basic waveform handling routines and structures for 
   libbpm

   The bpmwf sublibrary implements 3 waveform types doublewf_t, intwf_t and complexwf_t, 
   all of which are simple structure typedefs which hold the number of samples, the 
   sampling frequency and a pointer "wf" to the waveform. So the data array is accessible
   via doublewf_t::wf as a normal array of integers, doubles and complex_t 's. 

   @section wave_mem Memory management

   All have memory management routines (allocation/deletion) and routines to cast to 
   other times ( eg doublewf_t -> intwf_t or the other way around ). This can be done 
   either by filling existing waveforms ( convenient when you e.g. have already allocated
   memory and referenced it into a root branch ) or by having the casting routine allocate 
   memory itself and return a pointer to it. e.g:

   @code
   intwf_t *w = intwf_cast_new( doublewf_t *dw );
   @endcode

   this allocates memory for intwf_t and returns a pointer it, or

   @code
   intwf_cast( intwf_t *w, doublewf_t *dw );  
   @endcode

   this casts dw into existing intwf w.


   The sublibrary employs the sampling convention, where the sample is
   taken at the time index corresponding to 

   @code
   t = (double) i / sampling_freq
   @endcode

   @section wave_handle Waveform handling

   The sublibrary implements basic waveform hanlding like addition, subtraction,
    multiplication, division, biasing and scaling.

   Some advanced routines like differentiation, integration of the waveforms are also
   present. Also interpolation is impemented using various schemes which are more applicable 
   depending on the type of waveform : linear, parabolic : for non repeatative signals, 
   sinc and lanczos for repeatative signals (cfr. Shannon-Whittaker interpolation). 
   (thinking of cubic-spline as well... but not implemented yet). Using these interpolation 
   schemes, the sublibrary also implements resampling routines. 

   The complex waveforms have a set of routines to extract real/imag parts as well as
   phase and amplitude. Similar comments apply as for the casting routines, where the 
   "_new" versions allocate memory in the routine and return a pointer to it. 

   @section wave_fill Filling the waveforms

   The values of the waveforms can be set by either filling them from a given array of 
   values using e.g. 

   @code
   doublewf_setvalues( doublewf_t *w, double *a)
   @endcode

   or by calculating them from a function which returns the basic type of 
   the waveform.

   E.g. define a complex valued function in your code:

   @code
   complex_t csin( double t, int npars, double a ) {
     complex_t z
     // calculate a complex number z from the time t and parameters...
     return z;
   }
   @endcode

   which returns a complex value from the time \c t and having \c npar paramaters 
   \c a[0]  ... \c a[n-1]

   You can fill a waveform ( and so bascially sample the function at sampling frequency fs ) 
   by executing

   @code
   complexwf_setfunction( complexwf_t *z, &csin, npars, a )
   @endcode

   Also some routines are added to fill the waveforms with CW tones and decaying waves,
   along with some noise adding routines etc...


   @section wave_interpol Note on the interpolation options. 
   Here are some examples of the different interpolation options that one can 
   give to the doublewf/complexwf_getvalue() or _resample() routines.

   @image html linear_interpolation.gif "The WF_LINEAR interpolation option"
   @image html quadratic_interpolation.gif "The WF_QUADRATIC interpolation ption"
   @image html sinc_interpolation.gif "The WF_SINC interpolation option"
   @image html lanczos_interpolation.gif "The WF_LANCZOS interpolation option"

   @section wave_example For examples... 
   For examples on library use, please see the examples/wf 
   directory in the libbpm main tree...
   
   @section wave_todo Todo list
   - implement cubic spline interpolation ?
*/

/**
   @file
   @ingroup wave
   @brief   Simple waveform handling routines for libbpm
*/

/** @addtogroup wave */
/** @{ */

#ifndef BPMWAVE_H_LOADED__
#define BPMWAVE_H_LOADED__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "bpm/bpm_defs.h"
#include "bpm/bpm_units.h"
#include "bpm/bpm_messages.h"
#include "bpm/bpm_nr.h"

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */
#define WF_EPS          1.e-10 /**< A small number */
#define MAX_ALLOWED_NS  262144 /**< Maximum allowed number of samples (2^18)*/

#define WF_NEAREST    0x0000   /**< No interpolation, return nearest sample  */
#define WF_LINEAR     0x0001   /**< Perform linear interpolation in XXXwf_getsample() */
#define WF_QUADRATIC  0x0002   /**< Perform quadratic (parabolic) interpolation */
#define WF_SINC       0x0004   /**< signal reconstruction using sinc kernel (0..ns)*/     
#define WF_LANCZOS    0x0008   /**< signal reconstruction using lanczos kernel (a=3) */

/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

  /** Structure representing a waveform of doubles */
  typedef struct {
    int     ns; /**< The number of samples in the waveform */
    double  fs; /**< The sampling frequency */
    double *wf; /**< Pointer to an array of doubles which hold the samples */
  } doublewf_t;
  
  /** Structure representing a waveform of integers */
  typedef struct {
    int     ns; /**< The number of samples in the waveform */
    double  fs; /**< The sampling frequency */
    int    *wf; /**< Pointer to an array of integers which hold the samples */
  } intwf_t;
  
  /** Structure representing a waveform of complex numbers */
  typedef struct {
    int        ns; /**< The number of samples in the waveform */
    double     fs; /**< The sampling frequency */
    complex_t *wf; /**< Pointer to an array of integers which hold the samples */
  } complexwf_t;
  

  /** Structure with basic waveform statistics */
  typedef struct {
    int     imax;  /**< The sample nr of maximum of waveform */
    int     imin;  /**< The sample nr of minimum of waveform */
    double  max;   /**< The maximum value of waveform */
    double  min;   /**< The minimum value of waveform */
    double  mean;  /**< The mean of waveform  */
    double  rms;   /**< The rms of waveform */
  } wfstat_t;
  
/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */

  // wfstats_t handling routines -----------------------------------------------

  /** Reset the waveform statistics structure.
      @param s A pointer to a wfstat_t structure
      @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int wfstat_reset( wfstat_t *s );


  /** Prints the waveform statistics to the screen,
      @param of A filepointer
      @param s  A pointer to the waveform statistics structure
      @return void */
  EXTERN void wfstat_print( FILE *of, wfstat_t *s );
  

  // doublewf_t handling routines ----------------------------------------------

  /** Allocates memory for a new waveform of doubles 
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN doublewf_t* doublewf( int ns, double fs );

  
  /** Allocates memory for a new waveform of doubles and fills
      it with the sample time values
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN doublewf_t* doublewf_time_series( int ns, double fs );


  /** Allocates memory for a new waveform of doubles and fills
      it with sample numbers.
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN doublewf_t* doublewf_sample_series( int ns, double fs );


  /** Allocates memory for a new waveform of doubles and fills
      it with the frequency values
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN doublewf_t* doublewf_frequency_series( int ns, double fs );


  /** Fills the waveform of doubles with the values from the array x. 
      No check is performed whether x contains enough samples, the 
      user needs to be sure this is the case ! 
      @param w A pointer to the waveform of doubles
      @param x A pointer to the x values
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_setvalues( doublewf_t *w, double *x );

 
  /** Fills the waveform with values from the function wffun(), this function
      has to return a double from argument t ( time ) and has npars parameters
      given by the array *par. The function will be evaluated at the time t
      of each sample...
      @param w     A pointer to the waveform of doubles
      @param wffun A pointer to the function to fill the waveform with
      @param t     The time parameter in the function
      @param npars Number of parameters for the function
      @param par   Array of parameters for the function
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_setfunction( doublewf_t *w, 
				   double (*wffun)( double t, int, double* ),
				   int npars, double *par );
  

  /** Copies the values from existing waveform src into copy
      checks first whether the waveforms are compatible...
      This routine doesn't allocate memory internally and the
      waveforms should already have been created by the user...
      @param copy A pointer to the copy waveform
      @param src A pointer to the original waveform
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure.*/
  EXTERN int doublewf_copy( doublewf_t *copy, doublewf_t *src );


  /** Allocates memory and produces a copy of the waveform w;
      @param w A pointer to the original waveform
      @return A pointer to the copy of w */
  EXTERN doublewf_t* doublewf_copy_new( doublewf_t *w );


  /** Copies a subset from sample i1 to sample i2 ( inclusive ) to 
      the sub waveform from waveform w. The routine expects the 
      sub waveform to already exist with enough samples. ( this is
      not checked ! ) The sub->fs and sub->ns will be overwritten.
      @param sub Pointer to the waveform which will hold the subset
      @param w   Pointer to the original waveform
      @param i1  First sample of w to copy
      @param i2  Last sample of w to copy
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure.*/
  EXTERN int doublewf_subset( doublewf_t *sub, doublewf_t *w, int i1, int i2 );


  /** Resets the waveform of doubles to 0.
      @param w A pointer to the waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_reset( doublewf_t *w ); 


  /** Frees up the memory used by the waveform
      @param w A pointer to the waveform of doubles
      @return void */
  EXTERN void doublewf_delete( doublewf_t *w );
  

  /** Cast the waveform of doubles to a new waveform of integers.
      Memory is allocated inside this routine so the user just
      needs to have a inwf_t pointer ready. 
      @param w A pointer to the waveform of doubles
      @return A newly created intwf_t representation of the waveform
      of doubles  */
  EXTERN intwf_t* intwf_cast_new( doublewf_t *w );


  /** Cast the waveform of doubles to an already existing waveform 
      of integers.
      @param iw A pointer to an existing waveform of integers
      @param w  A pointer to the waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_cast( intwf_t *iw, doublewf_t *w );


  /** Checks compatiblity of the two waveforms, returns true 
      if the number of samples and the sampling frequencies match.
      For the sampling frequency, it is simply checked whether
      they match to WF_EPS.
      @param w1 A pointer to the first waveform of doubles
      @param w2 A pointer to the second waveform of doubles
      @return 1 if the waveforms match, 0 if not. */
  EXTERN int doublewf_compat( doublewf_t *w1, doublewf_t *w2 );

  
  /** Adds two waveforms of doubles w1+w2 sample per sample. 
      The result is  stored in w1.
      @param w1 A pointer to the first waveform of doubles
      @param w2 A pointer to the second waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_add( doublewf_t *w1, doublewf_t *w2 );


  /** Subtracts two waveforms of doubles w1-w2 sample per sample.
      The result is stored in w1.
      @param w1 A pointer to the first waveform of doubles
      @param w2 A pointer to the second waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_subtract( doublewf_t *w1, doublewf_t *w2 );


  /** Multiplies two waveforms of doubles w1*w2 sample per sample. 
      The result is stored in w1.
      @param w1 A pointer to the first waveform of doubles
      @param w2 A pointer to the second waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_multiply( doublewf_t *w1, doublewf_t *w2 );


  /** Divides two waveforms of doubles w1/w2 sample per sample. 
      The result is stored in w1. When w2[i] is 0, w1[i] will be
      set to 0. and a warning message is printed.
      @param w1 A pointer to the first waveform of doubles
      @param w2 A pointer to the second waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_divide( doublewf_t *w1, doublewf_t *w2 );


  /** Scales the waveform of doubles w by factor f.
      The result is stored in w. 
      @param f The scalefactor
      @param w A pointer to the waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_scale( double f, doublewf_t *w );


  /** Biases the waveform of doubles w by a constant c.
      The result is stored in w. 
      @param c The constant bias.
      @param w A pointer to the waveform of doubles
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_bias( double c, doublewf_t *w );


  /** Adds a cosine-like CW tone to the entire waveform. The sampling time
      is taken on the array index, so t=(double)i/w->fs.
      @param w     A pointer to the waveform structure
      @param amp   Amplitude of the CW tone
      @param phase Phase of the CW tone
      @param freq  Frequency of the CW tone
      @param phasenoise Sigma of the gaussian phasenoise 
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_add_cwtone( doublewf_t *w, double amp, double phase, double freq,
				  double phasenoise );

  
  /** Adds a decaying wave pulse to the waveform. The sampling time
      is taken on the array index, so t=(double)i/w->fs. The added signal
      is of the form :
      @f[ amp e^{- ( t - ttrig ) / tdcy } cos( 2 \pi freq ( t- ttrig ) + phase ) @f]
      If desired, phasenoise is added to the phase of the waveform. 
      @param w     A pointer to the waveform structure
      @param amp   Amplitude of the CW tone
      @param phase Phase of the CW tone
      @param freq  Frequency of the CW tone
      @param ttrig Trigger time of the pulse
      @param tdcy  Decay time of the pulse 
      @param phasenoise Sigma of the gaussian phasenoise 
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_add_dcywave( doublewf_t *w, double amp, double phase, double freq, 
				   double ttrig, double tdcy, double phasenoise );
  
  
  /** Adds gaussian amplitude noise to the waveform. 
      @param w     A pointer to the waveform structure
      @param sigma The gaussian sigma of the amplitude noise
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_add_ampnoise( doublewf_t *w, double sigma );

  

  /** Retrieves some basic statistics about the waveform of doubles in w,
      only considers samples between s0 and s1.
      @param w A pointer to the waveform structure
      @param s0 First sample to consider
      @param s1 Last sample to consider
      @param stats A filled wfstat_t structure is returned.
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_basic_stats( doublewf_t *w, int s0, int s1, wfstat_t *stats );


  /** Produce the derivative waveform for w : dw/dt.
      @param w A pointer to the waveform structure.
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int doublewf_derive( doublewf_t *w );

  /** Produce the integrated waveform for w : \int_0^t w(s)ds.
      @param w A pointer to the waveform structure.
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int doublewf_integrate( doublewf_t *w );


  /** Print the waveform to the filepointer
      @param of A filepointer, use stdout for the terminal
      @param w  A pointer to the waveform
      @return void
   */
  EXTERN void doublewf_print( FILE *of, doublewf_t *w );


  /** Return the value for the waveform at sample time t,
      according to the interpolation mode.
      @param w A pointer to the waveform structure
      @param t A time at which to sample the waveform
      @param mode Interpolation mode
      @return the value of the waveform at time t
   */
  EXTERN double doublewf_getvalue( doublewf_t *w, double t, unsigned int mode );
  

  /** Resamples the waveform w1 into w2 with new fs sampling frequency
      This routine recalculates the correct number of samples required.
      However the user needs to make sure that there are enough samples
      in w2 available as this is not checked. The w2->ns value will be
      overwritten with the correct amount. 
      The routine checkes whether the maximum allowed number of samples
      is not exceeded to avoid memory problems.
      @param w A pointer to the waveform structure
      @param t A time at which to sample the waveform
      @param mode Interpolation mode
      @return the value of the waveform at time t
   */
  EXTERN int doublewf_resample( doublewf_t *w2, double fs, doublewf_t *w1, unsigned int mode );



  // intwf_t handling routines ----------------------------------------------


  /** Allocates memory for a new waveform of integers
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN intwf_t* intwf( int ns, double fs );
  
  /** Allocates memory for a new waveform of integers and fills
      it with sample numbers.
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN intwf_t* intwf_sample_series( int ns, double fs );


  /** Fills the waveform of integers with the values from the array x. 
      No check is performed whether x contains enough samples, the 
      user needs to be sure this is the case ! 
      @param w A pointer to the waveform of integers
      @param x A pointer to the x values
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_setvalues( intwf_t *w, int *x );

 
  /** Fills the waveform with values from the function wffun(), this function
      has to return a double from argument t ( time ) and has npars parameters
      given by the array *par. The function will be evaluated at the time t
      of each sample...
      @param w     A pointer to the waveform of integers
      @param wffun A pointer to the function to fill the waveform with
      @param t     The time parameter in the function
      @param npars Number of parameters for the function
      @param par   Array of parameters for the function
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_setfunction( intwf_t *w, 
				int (*wffun)( double t, int, double* ),
				int npars, double *par );
  

  /** Copies the values from existing waveform src into copy
      checks first whether the waveforms are compatible...
      This routine doesn't allocate memory internally and the
      waveforms should already have been created by the user...
      @param copy A pointer to the copy waveform
      @param src A pointer to the original waveform
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure.*/
  EXTERN int intwf_copy( intwf_t *copy, intwf_t *src );


  /** Allocates memory and produces a copy of the waveform w;
      @param w A pointer to the original waveform
      @return A pointer to the copy of w */
  EXTERN intwf_t* intwf_copy_new( intwf_t *w );


  /** Copies a subset from sample i1 to sample i2 ( inclusive ) to 
      the sub waveform from waveform w. The routine expects the 
      sub waveform to already exist with enough samples. ( this is
      not checked ! ) The sub->fs and sub->ns will be overwritten.
      @param sub Pointer to the waveform which will hold the subset
      @param w   Pointer to the original waveform
      @param i1  First sample of w to copy
      @param i2  Last sample of w to copy
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure.*/
  EXTERN int intwf_subset( intwf_t *sub, intwf_t *w, int i1, int i2 );


  /** Resets the waveform of integers to 0.
      @param w A pointer to the waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_reset( intwf_t *w ); 


  /** Frees up the memory used by the waveform
      @param w A pointer to the waveform of integers
      @return void */
  EXTERN void intwf_delete( intwf_t *w );
  

  /** Cast the waveform of integers to a new waveform of doubles.
      Memory is allocated inside this routine so the user just
      needs to have a inwf_t pointer ready. 
      @param w A pointer to the waveform of integers
      @return A newly created doublewf_t representation of the waveform
      of integers  */
  EXTERN doublewf_t* doublewf_cast_new( intwf_t *w );


  /** Cast the waveform of integers to an already existing waveform 
      of doubles.
      @param iw A pointer to an existing waveform of integers
      @param w  A pointer to the waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int doublewf_cast( doublewf_t *w, intwf_t *iw );


  /** Checks compatiblity of the two waveforms, returns true 
      if the number of samples and the sampling frequencies match.
      For the sampling frequency, it is simply checked whether
      they match to WF_EPS.
      @param w1 A pointer to the first waveform of integers
      @param w2 A pointer to the second waveform of integers
      @return 1 if the waveforms match, 0 if not. */
  EXTERN int intwf_compat( intwf_t *w1, intwf_t *w2 );

  
  /** Adds two waveforms of integers w1+w2 sample per sample. 
      The result is  stored in w1.
      @param w1 A pointer to the first waveform of integers
      @param w2 A pointer to the second waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_add( intwf_t *w1, intwf_t *w2 );


  /** Subtracts two waveforms of integers w1-w2 sample per sample.
      The result is stored in w1.
      @param w1 A pointer to the first waveform of integers
      @param w2 A pointer to the second waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_subtract( intwf_t *w1, intwf_t *w2 );


  /** Multiplies two waveforms of integers w1*w2 sample per sample. 
      The result is stored in w1.
      @param w1 A pointer to the first waveform of integers
      @param w2 A pointer to the second waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_multiply( intwf_t *w1, intwf_t *w2 );


  /** Divides two waveforms of integers w1/w2 sample per sample. 
      The result is stored in w1. When w2[i] is 0, w1[i] will be
      set to 0. and a warning message is printed.
      @param w1 A pointer to the first waveform of integers
      @param w2 A pointer to the second waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_divide( intwf_t *w1, intwf_t *w2 );


  /** Scales the waveform of integers w by factor f.
      The result is stored in w. 
      @param f The scalefactor
      @param w A pointer to the waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_scale( int f, intwf_t *w );


  /** Biases the waveform of integers w by a constant c.
      The result is stored in w. 
      @param c The constant bias.
      @param w A pointer to the waveform of integers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_bias( int c, intwf_t *w );


  /** Adds a cosine-like CW tone to the entire waveform. The sampling time
      is taken on the array index, so t=(double)i/w->fs.
      @param w     A pointer to the waveform structure
      @param amp   Amplitude of the CW tone
      @param phase Phase of the CW tone
      @param freq  Frequency of the CW tone
      @param phasenoise Sigma of the gaussian phasenoise 
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_add_cwtone( intwf_t *w, double amp, double phase, double freq,
			       double phasenoise );

  
  /** Adds a decaying wave pulse to the waveform. The sampling time
      is taken on the array index, so t=(double)i/w->fs. The added signal
      is of the form :
      @f[ amp e^{- ( t - ttrig ) / tdcy } cos( 2 \pi freq ( t- ttrig ) + phase ) @f]
      If desired, phasenoise is added to the phase of the waveform. 
      @param w     A pointer to the waveform structure
      @param amp   Amplitude of the CW tone
      @param phase Phase of the CW tone
      @param freq  Frequency of the CW tone
      @param ttrig Trigger time of the pulse
      @param tdcy  Decay time of the pulse 
      @param phasenoise Sigma of the gaussian phasenoise 
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_add_dcywave( intwf_t *w, double amp, double phase, double freq, 
				double ttrig, double tdcy, double phasenoise );
  
  
  /** Adds gaussian amplitude noise to the waveform. 
      @param w     A pointer to the waveform structure
      @param sigma The gaussian sigma of the amplitude noise
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_add_ampnoise( intwf_t *w, double sigma );

  

  /** Retrieves some basic statistics about the waveform of integers in w,
      only considers samples between s0 and s1.
      @param w A pointer to the waveform structure
      @param s0 First sample to consider
      @param s1 Last sample to consider
      @param stats A filled wfstat_t structure is returned.
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int intwf_basic_stats( intwf_t *w, int s0, int s1, wfstat_t *stats );


  /** Produce the derivative waveform for w : dw/dt.
      @param w A pointer to the waveform structure.
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int intwf_derive( intwf_t *w );

  /** Produce the integrated waveform for w : \int_0^t w(s)ds.
      @param w A pointer to the waveform structure.
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure */
  EXTERN int intwf_integrate( intwf_t *w );


  /** Print the waveform to the filepointer
      @param of A filepointer, use stdout for the terminal
      @param w  A pointer to the waveform
      @return void
   */
  EXTERN void intwf_print( FILE *of, intwf_t *w );


  /** Return the value for the waveform at sample time t,
      according to the interpolation mode.
      @param w A pointer to the waveform structure
      @param t A time at which to sample the waveform
      @param mode Interpolation mode
      @return the value of the waveform at time t
   */
  EXTERN int intwf_getvalue( intwf_t *w, double t, unsigned int mode );
  

  /** Resamples the waveform w1 into w2 with new fs sampling frequency
      This routine recalculates the correct number of samples required.
      However the user needs to make sure that there are enough samples
      in w2 available as this is not checked. The w2->ns value will be
      overwritten with the correct amount. 
      The routine checkes whether the maximum allowed number of samples
      is not exceeded to avoid memory problems.
      @param w A pointer to the waveform structure
      @param t A time at which to sample the waveform
      @param mode Interpolation mode
      @return the value of the waveform at time t
   */
  EXTERN int intwf_resample( intwf_t *w2, double fs, intwf_t *w1, unsigned int mode );



  // complexwf_t handling routines ----------------------------------------------

  /** Allocates memory for a new waveform of complex numbers
      @param ns The number of samples in the waveform
      @param fs The sampling frequency of the waveform
      @return A pointer to the allocated waveform structure */
  EXTERN complexwf_t* complexwf( int ns, double fs );

  /** Allocates memory and produces a copy of the complex waveform w;
      @param w A pointer to the original waveform
      @return A pointer to the copy of w */
  EXTERN complexwf_t* complexwf_copy_new( complexwf_t *w );

  /** Copies the values from existing complex waveform src into copy
      checks first whether the waveforms are compatible...
      This routine doesn't allocate memory internally and the
      waveforms should already have been created by the user...
      @param copy A pointer to the copy waveform
      @param src A pointer to the original waveform
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure.*/
  EXTERN int complexwf_copy( complexwf_t *copy, complexwf_t *src );

 /** Copies a subset from sample i1 to sample i2 ( inclusive ) to 
      the sub waveform from complex waveform w. The routine expects the 
      sub waveform to already exist with enough samples. ( this is
      not checked ! ) The sub->fs and sub->ns will be overwritten.
      @param sub Pointer to the waveform which will hold the subset
      @param w   Pointer to the original waveform
      @param i1  First sample of w to copy
      @param i2  Last sample of w to copy
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure.*/
  EXTERN int complexwf_subset( complexwf_t *sub, complexwf_t *w, int i1, int i2 );


  /** Fills the complex waveform with the values from the array x. 
      No check is performed whether x contains enough samples, the 
      user needs to be sure this is the case ! 
      @param w A pointer to the waveform of complex numbers
      @param x A pointer to the complex x values
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_setvalues( complexwf_t *w, complex_t *x );

  /** Fills the waveform with values from the function wffun(), this function
      has to return a complex_t from argument t ( time ) and has npars parameters
      given by the array *par. The function will be evaluated at the time t
      of each sample...
      @param w     A pointer to the waveform of complex numbers
      @param wffun A pointer to the function to fill the waveform with
      @param t     The time parameter in the function
      @param npars Number of parameters for the function
      @param par   Array of parameters for the function
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_setfunction( complexwf_t *w, 
				    complex_t (*wffun)( double, int, double* ),
				    int npars, double *par );

  /** Resets the waveform of complex numbers to 0+0i
      @param w A pointer to the complex waveform
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_reset( complexwf_t *w );

  /** Frees up the memory used by the waveform
      @param w A pointer to the waveform of complex numbers
      @return void */
  EXTERN void complexwf_delete( complexwf_t *w );

  /** Checks compatiblity of the two waveforms, returns true 
      if the number of samples and the sampling frequencies match.
      For the sampling frequency, it is simply checked whether
      they match to WF_EPS.
      @param w1 A pointer to the first waveform of complex numbers
      @param w2 A pointer to the second waveform of complex numbers
      @return 1 if the waveforms match, 0 if not. */
  EXTERN int complexwf_compat( complexwf_t *w1, complexwf_t *w2 );

  /** Adds two waveforms of complex numbers w1+w2 sample per sample. 
      The result is  stored in w1.
      @param w1 A pointer to the first waveform of complex numbers
      @param w2 A pointer to the second waveform of comlex numbers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_add( complexwf_t *w1, complexwf_t *w2 );

  /** Subtracts two waveforms of complex numbers w1-w2 sample per sample. 
      The result is  stored in w1.
      @param w1 A pointer to the first waveform of complex numbers
      @param w2 A pointer to the second waveform of comlex numbers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_subtract( complexwf_t *w1, complexwf_t *w2 );

  /** Multiplies two waveforms of complex numbers w1*w2 sample per sample. 
      The result is  stored in w1.
      @param w1 A pointer to the first waveform of complex numbers
      @param w2 A pointer to the second waveform of comlex numbers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_multiply( complexwf_t *w1, complexwf_t *w2 );

  
  /** Divides two waveforms of complex numbers w1/w2 sample per sample. 
      The result is  stored in w1.
      @param w1 A pointer to the first waveform of complex numbers
      @param w2 A pointer to the second waveform of comlex numbers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_divide( complexwf_t *w1, complexwf_t *w2 );


  /** Scales the waveform of complex numbers w with complex factor f
      The result is  stored in w.
      @param f  The complex scaling factor
      @param w  A pointer to the waveform of comlex numbers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_scale( complex_t f, complexwf_t *w );


  /** Biases the waveform of complex numbers w with complex constant c
      The result is  stored in w.
      @param c  The complex constant
      @param w  A pointer to the waveform of comlex numbers
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_bias( complex_t c, complexwf_t *w );


  /** Adds a CW tone to the entire waveform. The sampling time
      is taken on the array index, so t=(double)i/w->fs.
      The real part will have the cos-like waveform, the imaginary part
      the sin-like waveform.
      @param w     A pointer to the complex waveform structure
      @param amp   Amplitude of the CW tone
      @param phase Phase of the CW tone
      @param freq  Frequency of the CW tone
      @param phasenoise Sigma of the gaussian phasenoise 
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_add_cwtone( complexwf_t *w, double amp, double phase, double freq,
				   double phasenoise );
  

  /** Adds a decaying wave pulse to the waveform. The sampling time
      is taken on the array index, so t=(double)i/w->fs. The added signal
      is of the form :
      @f[ amp e^{- ( t - ttrig ) / tdcy } sin( 2 \pi freq ( t- ttrig ) + phase ) @f]
      The real part will have the cos-like component, the imaginary part
      the sin-like component. If desired, phasenoise is added to the phase of the 
      waveform. 
      @param w     A pointer to the waveform structure
      @param amp   Amplitude of the CW tone
      @param phase Phase of the CW tone
      @param freq  Frequency of the CW tone
      @param ttrig Trigger time of the pulse
      @param tdcy  Decay time of the pulse 
      @param phasenoise Sigma of the gaussian phasenoise 
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */

  EXTERN int complexwf_add_dcywave( complexwf_t *w, double amp, double phase, double freq, 
				    double ttrig, double tdcy, double phasenoise );
 
  /** Adds uncorrelated gaussian amplitude noise with uniformly distributed random 
      phase to the complex the waveform. 
      @param w     A pointer to the complex waveform structure
      @param sigma The gaussian sigma of the amplitude noise, phase is uniform over 2pi
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_add_noise( complexwf_t *w, double sigma );


  /** Adds pure gaussian amplitude noise to the complex waveform and leaves 
      the phase untouched
      @param w     A pointer to the complex waveform structure
      @param sigma The gaussian sigma of the amplitude noise
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_add_ampnoise( complexwf_t *w, double sigma );


  /** Adds pure gaussian phase noise to the complex waveform and leaves 
      the amplitude untouched
      @param w     A pointer to the complex waveform structure
      @param sigma The gaussian sigma of the phase noise
      @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure. */
  EXTERN int complexwf_add_phasenoise( complexwf_t *w, double sigma );


  /** Print the waveform to the filepointer
      @param of A filepointer, use stdout for the terminal
      @param w  A pointer to the waveform
      @return void
  */
  EXTERN void complexwf_print( FILE *of, complexwf_t *w );

  /**
     Gets the real part of the comlex waveform into the waveform of doubles.
     The doublewf needs to be allocated by the user beforehand and have the 
     same number of samples as the complex waveform.
     @param re A pointer to the waveform of doubles which will store the real part
     @param z  A pointer to the complex waveform
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int complexwf_getreal( doublewf_t *re, complexwf_t *z );

  /**
     Gets the imaginary part of the comlex waveform into the waveform of doubles.
     The doublewf needs to be allocated by the user beforehand and have the 
     same number of samples as the complex waveform.
     @param im A pointer to the waveform of doubles which will store the imaginary part
     @param z  A pointer to the complex waveform
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int complexwf_getimag( doublewf_t *im, complexwf_t *z );

  /**
     Gets the amplitude of the comlex waveform into the waveform of doubles.
     The doublewf needs to be allocated by the user beforehand and have the 
     same number of samples as the complex waveform.
     @param im A pointer to the waveform of doubles which will store the amplitude
     @param z  A pointer to the complex waveform
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int complexwf_getamp( doublewf_t *r, complexwf_t *z );

  /**
     Gets the phase of the comlex waveform into the waveform of doubles.
     The doublewf needs to be allocated by the user beforehand and have the 
     same number of samples as the complex waveform. The phase is normalised
     between [0,2pi[.
     @param im A pointer to the waveform of doubles which will store the phase
     @param z  A pointer to the complex waveform
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int complexwf_getphase( doublewf_t *theta, complexwf_t *z );

  /** Retrieves the real part of the complex waveform in a newly allocated
      waveform of doubles. Memory on the heap is allocated inside this routine, 
      the user has to deal with deal with freeing it him/her self.
      @param z A pointer to the complex waveform
      @return A pointer to the allocated waveform of doubles containing the real
      part of z. */
  EXTERN doublewf_t* complexwf_getreal_new( complexwf_t *z );

  /** Retrieves the imaginary part of the complex waveform in a newly allocated
      waveform of doubles. Memory on the heap is allocated inside this routine, 
      the user has to deal with deal with freeing it him/her self.
      @param z A pointer to the complex waveform
      @return A pointer to the allocated waveform of doubles containing the imaginary
      part of z. */
  EXTERN doublewf_t* complexwf_getimag_new( complexwf_t *z );

  /** Retrieves the amplitude of the complex waveform in a newly allocated
      waveform of doubles. Memory on the heap is allocated inside this routine, 
      the user has to deal with deal with freeing it him/her self.
      @param z A pointer to the complex waveform
      @return A pointer to the allocated waveform of doubles containing the amplitude
      of z. */
  EXTERN doublewf_t* complexwf_getamp_new( complexwf_t *z );

  /** Retrieves the phase of the complex waveform in a newly allocated
      waveform of doubles. Memory on the heap is allocated inside this routine, 
      the user has to deal with deal with freeing it him/her self. The phase is
      normalised between [0,2pi[.
      @param z A pointer to the complex waveform
      @return A pointer to the allocated waveform of doubles containing the phase
      of z. */
  EXTERN doublewf_t* complexwf_getphase_new( complexwf_t *z );

  /**
     Set the real part of the complex waveform z to re. The complexwf needs 
     to be allocated by the user beforehand and have the same number of samples 
     as the double waveform.
     @param z  A pointer to the complex waveform
     @param re A pointer to a waveform of double containing the real part
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int complexwf_setreal( complexwf_t *z, doublewf_t *re );

  /**
     Set the imaginary part of the complex waveform z to im. The complexwf needs 
     to be allocated by the user beforehand and have the same number of samples 
     as the double waveform.
     @param z  A pointer to the complex waveform
     @param re A pointer to a waveform of double containing the imaginary part
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure */
  EXTERN int complexwf_setimag( complexwf_t *z, doublewf_t *im );


  // Sample / Frequency / Time conversionn routines -----------------------------
  EXTERN int time_to_sample( double fs, int ns, double t, int *iS );

  EXTERN int freq_to_sample( double fs, int ns, double f, int *iS );

  EXTERN int sample_to_time( double fs, int ns, int iS, double *t );

  EXTERN int sample_to_freq( double fs, int ns, int iS, double *f );

#ifdef __cplusplus
}
#endif

#endif /* #ifndef __BPMWAVE_H_LOADED */
/** @} */
/* ========================== end of file 'bpm_wave.h' ==================== */
