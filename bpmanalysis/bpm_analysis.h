/**
   @file    
   @ingroup analysis Analysis routines
   @brief   libbpm analysis routines

   This header contains definitions for the libbpm BPM data analysis routines.
   These mainly are the SVD and resolution/residual calculation routines along
   with the definition of an analysis cut function...
*/

/** @addtogroup analysis */
/** @{ */

#ifndef BPMANALYSIS_H__
#define BPMANALYSIS_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <math.h>
#include <bpm/bpm_defs.h>
#include <bpm/bpm_interface.h>

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */

#define BPM_GOOD_EVENT    1   /**< A good event */
#define BPM_BAD_EVENT     0   /**< A bad event */

#define ANA_SVD_TILT      0x00000001   /**< Include tilts in the SVD */
#define ANA_SVD_NOTILT    0x00000000   /**< Don't include tilts in the SVD */

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
     Set the cut function

     @param cutfn a pointer to the cut function with a bpmproc_t as argument
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */

  EXTERN int ana_set_cutfn( int (*cutfn)( bpmproc_t *proc ) );


  /**
     Perform the SVD on the given data and return the coefficients. The index 0 bpmconf
     is the bpm to be regressed against and the remainder are put into the regression.
     The coeffs array must be valid up to the number of arguments approapriate to mode.
     
     @param proc      pointer to the the processed bpm databuffer
     @param num_bpms  the number of bpms in the array
     @param num_svd   number of svd constants
     @param total_num_evts total number of events in the buffer
     @param coeffs    the array of correlation coefficients that is returned
     @param mode      mode option: take tilts into account in the SVD ?
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int ana_get_svd_coeffs( bpmproc_t **proc, int num_bpms, int num_svd, 
				 int total_num_evts, double *coeffs, int mode );
  
  
  /**
     Calculate the mean and rms of the residual fomr the given events. Note that the
     mode and svd coefficients must 'match' as with ana_get_svd_coeffs()

     @param proc      pointer to the the processed bpm databuffer
     @param num_bpms  the number of bpms in the array
     @param num_evts  total number of events in the buffer
     @param coeffs    the array of correlation coefficients
     @param mode      mode option: take tilts into account in the SVD ?
     @param mean      the returned mean
     @param rms       the returned rms
  */
  EXTERN int ana_compute_residual( bpmproc_t **proc, int num_bpms, int num_evts,
				   double *coeffs, int mode, double *mean, double *rms );

  /**
     The default cut function if people cut be bothered to do their own :)

     @param  proc the event to decide
     @return BPM_GOOD_EVENT if the event is good, BPM_BAD_EVENT if it isn't
  */           
  EXTERN int ana_def_cutfn( bpmproc_t *proc );


  /**
     A user cut function to allow cuts to be applied while selecting events for SVD, etc.
  */
  EXTERN int (*ana_cutfn)( bpmproc_t *proc );


#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMANALYSIS_H__ */
/*@}*/
/* ================================ end of file ============================= */
