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
   @file
   @ingroup rf
   @brief   libbpm rf simulation routines

   The header file for RF routines

   Need to check in how far these routines are redundant, bpmdsp can 
   replace most of the filtering routines here !
*/

/** @addtogroup rf */
/** @{ */

#ifndef BPMRF_H__
#define BPMRF_H__

/* -----------------------------------------------------------------------------
   includes
   -------------------------------------------------------------------------- */
#include <math.h>
#include <bpm/bpm_defs.h>
#include <bpm/bpm_interface.h>
#include <bpm/bpm_wf.h>

/* -----------------------------------------------------------------------------
   macro definitions
   -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
   typedefs, enums and other declarations
   -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------------
   function prototypes and declarations
   -------------------------------------------------------------------------- */
  /**
     Numer of samples in the rf waveform representations, default value 
     is 2^16 = 65536
   */
  EXTERN int      rf_nsamples; 

  /**
     Effective sampling frequency for the rf waveform representations, 
     default value is 20 GHz
  */
  EXTERN double   rf_samplefreq; 

  EXTERN int      rf_setup( int nsamples, double sfreq );

  /**
     Rectifies the given waveform assuming a single diode
     @param D the rectified signal
     @param RF the complex waveform to rectify
     @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int rf_rectify( doublewf_t *D, complexwf_t *RF );
  
  /**
    Generates an LO waveform
    @param amp amplitude of the LO signal in Volts
    @param lofreq LO frequency
    @type locked or freerunning oscillator
    @phase phase of the signal, ignored if type is not "locked"
    @phasenoise phase noise to be added to the waveform
    @param LO generated waveform
    @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int rf_addLO( double amp, double lofreq, enum bpmphase_t type,
		       double phase, double phasenoise, doublewf_t *LO );

  /**
    Simulates an ideal mixer
    @param RF signal to mix
    @param LO local oscillator signal to mix with
    @param IF resulting signal containing the up and down
    converted terms
    @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int rf_mixer( doublewf_t *RF_Re, doublewf_t *LO, doublewf_t *IF );

  /**
     Amplifies the signal by the level dB. The voltage gain is calculated:
     @f[ gain = \sqrt{ 10^{ \frac{db}{20} } } @f]
     @param RF waveform to be processed
     @param dB gain (or attenuation) in dB
     @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int rf_amplify( doublewf_t *RF, double dB );

  /**
     Amplifies the signal by the level dB. The voltage gain is calculated:
     @f[ gain = \sqrt{ 10^{ \frac{db}{20} } } @f]
     @param RF waveform to be processed
     @param dB gain (or attenuation) in dB
     @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int rf_amplify_complex( complexwf_t *RF, double dB );

  /**
     Rotates the phase of the signal by the amount specified
     @param RF waveform to be processed
     @param rotation phase rotation in degrees
     @returns BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int rf_phase_shifter( complexwf_t *RF, double rotation );

#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMRF_H__ */
/** @} */
/* ================================ end of file ============================= */
