/** 
    @file
    @ingroup calib
    @brief   calibration routines

    This header contains some BPM calibration routines
*/

/** @addtogroup calib */
/** @{ */

#ifndef BPMCALIBRATION_H__
#define BPMCALIBRATION_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <math.h>
#include <bpm/bpm_defs.h>
#include <bpm/bpm_interface.h>

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */

  /**
     This routine basically defines the calibration steps and returns them into
     the array of beam structures. It needs an array of processed waveform structures,
     of dimension npulses from a single BPM. From this it determines the corresponding 
     corrector/mover steps and puts them back into the array of beam structures given 
     the bpm configurations.  

     Startpulse and stoppulse have to be in the first and last calib steps &
     will need some extensive error checking for e.g. missed calibration steps...

     NOTE: This is not definitive yet - more checking, etc. required! 
            - DDC or FIT?
            - Sign errors? 
            - not robust to missing steps

     @param proc       array of processed waveforms for a single bpm, so array of pulses
     @param cnf        array of bpm configuration structures
     @param npulses    number of pulses in the calibration
     @param startpulse start of calibration range
     @param stoppulse  stop of calibration range
     @param angle
     @param startpos   start position of calibration
     @param endpos     end position of calibration
     @param num_steps  number of calibration steps
     @param bunch      the returned bunchconf array which represents where the beam is 
                       supposed to be in each bpm during each calibration step

     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int setup_calibration( bpmconf_t *cnf, bpmproc_t *proc, int npulses, int startpulse, 
				int stoppulse, double angle, double startpos, double endpos, 
				int num_steps, bunchconf_t *bunch );

  
  /**
     Gets the calibration constants from an array of npulses of beam positions
     and processed waveform structures and returns an updated calibration structure.
     Note that this routine updates the IQ phase, the position scale and the tilt scale 
     but DOES NOT touch the frequency, decay time or the t0Offset.

     @param  bpm     Bpm structures
     @param  bunch   An array of bunch structures, one for each pulse, so essentially 
                     this corresponds to where we expect the beam to be in each pulse,
                     so representing corrector positions or mover positions. This information 
                     should be filled by the routine setup_calibration( ... )
     @param  proc    An array of processed waveforms, one for each pulse, which correspond
                     to calculated positions that were calculated using IQ phase = 0 and
                     scales equal to 1.
     @param  npulses The number of pulses in the arrays 
     @param *cal     The returned calibration structure for the BPM that was calibrated 
     @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure 
  */ 
  EXTERN int calibrate( bpmconf_t *bpm, bunchconf_t *bunch, bpmproc_t *proc, int npulses, 
			bpmcalib_t *cal );


#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMCALIBRATION_H__ */
/*@}*/
/* ================================ end of file ============================= */
