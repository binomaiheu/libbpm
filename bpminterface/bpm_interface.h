/**
   @file    
   @ingroup interface
   @brief    Front end interface structure definitions and handlers.

   This header contains the front-end interface structures and handlers
   for libbpm. They define a set of user friendly structures like
   bpmconf_t, bpmcalib_t, beamconf_t etc... to work with the bpm data.
*/

/** @addtogroup interface */
/** @{ */

#ifndef BPMINTERFACE_H__
#define BPMINTERFACE_H__

/* ---------------------------------------------------------------------------------------
// includes
// ------------------------------------------------------------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bpm/bpm_defs.h>
#include <bpm/bpm_wf.h>
#include <bpm/bpm_dsp.h>

/* ----------------------------------------------------------------------------------------
// macro definitions
// ------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------
// typedefs, enums and other declarations
// ------------------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

  /**
     BPM cavity ( of better signal ) type
  */
  enum bpmtype_t  { 
    diode,           /**< rectified bpm signal ( trigger pulse )                         */
    monopole,        /**< reference cavity signal ( monopole )                           */
    dipole           /**< position sentivive cavity signal ( dipole )                    */
  };

  /**
     Diode behavior type
  */
  enum triggertype {
    positive,        /**< Positive half-period of the waveform is detected               */
    negative,        /**< Negative half-period of the waveform is detected               */
    bipolar          /**< The both half-periods are detected                             */
  };
 
  /**
     BPM polarisation plane, basically a difficult way to say x or y ;)
  */
  enum bpmpol_t   { 
    horiz,           /**< Horizontal plane, or x in most cases                           */
    vert             /**< Vertical plane, or y in most cases                             */
  };
  
  /**
     BPM electronics phase lock type
  */
  enum bpmphase_t { 
    randomised,      /**< unlocked phase                                                 */
    locked           /**< locked phase                                                   */
  };
  

  typedef struct bpmconf   bpmconf_t;    /**< type definition for BPM configuration      */
  typedef struct bpmcalib  bpmcalib_t;   /**< type definition for calibrations           */
  typedef struct bpmproc   bpmproc_t;    /**< type definition for processed BPM signals  */
  typedef struct beamconf  beamconf_t;   /**< type definition for beam configurations    */
  typedef struct bunchconf bunchconf_t;  /**< type definition for bunch configurations   */
  typedef struct bpmmode   bpmmode_t;   
  typedef struct rfmodel   rfmodel_t;

  typedef enum   triggertype triggertype_t; 

  /**
     Structure containing the BPM configuration
  */
  struct bpmconf {
    char     name[20];              /**< a BPM should have a name                        */
    
    enum bpmtype_t  cav_type;         /**< BPM type                                      */
    enum bpmpol_t   cav_polarisation; /**< BPM polarisation                              */
    enum bpmphase_t cav_phasetype;    /**< BPM phase type                                */

    rfmodel_t *cav_model;
    double     cav_length;          /**< length of the cavity                            */
    double     cav_freq;            /**< cavity freq (MHz)                               */
    double     cav_decaytime;       /**< cavity decay time (microsec)                    */
    double     cav_phase;           /**< phase advance wrt. reference (fixed or random)  */
    double     cav_iqrotation;      /**< cavity IQ rotation                              */
    double     cav_chargesens;      /**< charge sensitivity (volt/nC)                    */
    double     cav_possens;         /**< pos sensitivity at 1.6nC charge (volt/micron)   */
    double     cav_tiltsens;        /**< tilt sensitivity at 1.6nC charge (volt/micron)  */
    
    double     rf_LOfreq;           /**< LO frequency to mix down with (in MHz)          */
    
    // digitisation configuration parameters
    double     digi_trigtimeoffset; /**< time (usec) to offset bunch arrival times by    */
    double     digi_freq;           /**< digitization frequency (MHz)                    */
    int        digi_nbits;          /**< number of bits in ADC for digitisation          */
    int        digi_nsamples;       /**< number of samples in ADC digitisation           */
    double     digi_ampnoise;       /**< amplitude noise in ADC channels (pedestal width)*/
    int        digi_voltageoffset;  /**< voltage offset (pedestal position) in counts    */
    double     digi_phasenoise;     /**< phase noise                                     */
    

    /* ------------------- normal waveform processing parameters ----------------------- */
    double     t0;                  /**< start time of pulse                             */

    // ddc processing configuration parameters
    double     ddc_freq;            /**< Frequency of downmixed waveform (MHz)           */
    double     ddc_tdecay;          /**< Decay time (usec)                               */
    double     ddc_tOffset;         /**< Always have offset from t0 for sampling !!!     */
    filter_t  *ddc_filter;          /**< DDC 2 omega filter                              */
    
    // fit processing configuration parameters
    double     fit_inifreq;         /**< Initial frequency for fitting                   */
    double     fit_initdecay;       /**< Initial decay time for fitting                  */
    double     fit_tOffset;         /**< Offset from t0 to start fitting                 */

    /* -------------------- caltone waveform processing parameters --------------------- */
    double     ddc_ct_freq;        /**< caltone frequency for the ddc algorithm          */
    filter_t  *ddc_ct_filter;      /**< filter for the caltone ddc                       */
    int        ddc_ct_iSample;     /**< sample number to sample from ddc for amp/phase   */


    /* -------------------- geometry of the BPM in the beamline ------------------------ */
    double     geom_pos[3];         /**< position of the BPM in the beamline             */
    double     geom_tilt[3];        /**< tilt of the BPM (0: xrot, 1: yrot, 2: zrot)     */

    /* --------------------- indices to use for user programs -------------------------- */
    int        ref_idx;             /**< reference cavity index for this BPM             */
    int        diode_idx;           /**< reference diode index for this BPM              */
    int        forced_trigger;      /**< this cavity is abused as trigger signal         */

    /* ------------------------- some processing buffers ------------------------------- */
    doublewf_t* ddc_buffer_re;      /**< pointer to a doublewf_t buffer                  */
    doublewf_t* ddc_buffer_im;      /**< pointer to a doublewf_t buffer                  */
  };

  /**
     A structure containing the calibration information : purely calibration !
  */
  struct bpmcalib {
    // ddc calibration constants
    double    ddc_IQphase;     /**< processed IQ phase for the ddc routine               */
    double    ddc_posscale;    /**< processed position scale for the ddc routine         */
    double    ddc_slopescale;  /**< processed slope scale for the fit routine            */
    double    ddc_ct_amp;      /**< calibration tone amplitude at time of calibration    */
    double    ddc_ct_phase;    /**< calibration tone phase at time of calibration        */

    // fit calibration constants
    double    fit_IQphase;     /**< processed IQ phase for the fit routine               */
    double    fit_posscale;    /**< position scale for the fit routine                   */
    double    fit_slopescale;  /**< slope scale for the fit routine                      */
    double    fit_ct_amp;      /**< calibration tone amplitude at time of calibration    */
    double    fit_ct_phase;    /**< calibration tone phase at time of calibration        */
  };

  /**
     A structure containing the processed waveform information
  */
  struct bpmproc {
    double   ampnoise;        /**< calculated (processed) amplitude noise                */
    double   voltageoffset;   /**< calculated voltage offset                             */
    
    double   t0;              /**< trigger t0 for int, copied from bpmconf_t::t0 for ext */
    
    int      saturated;       /**< this signal was saturated                             */
    int      iunsat;          /**< the last unsaturated sample index                     */

    complexwf_t *dc;          /**< The signal's DC waveform                              */
    complexwf_t *ft;          /**< The signal's fourier transform                        */
    
    // FFT parameters
    int      fft_success;     /**< do we have proper fft info ?                          */
    double   fft_amp;         /**< amplitude of fft                                      */
    double   fft_freq;        /**< frequency obtained from fft (MHz)                     */
    double   fft_tdecay;      /**< decay time obtained from fft (usec)                   */
    double   fft_offset;      /**< offset of fft in fit                                  */
    
    // DDC parameters
    int      ddc_success;     /**< do we have proper ddc info ?                          */
    double   ddc_tSample;     /**< time at which the ddc was sampled, t0+t0Offset        */
    int      ddc_iSample;     /**< index of sample at which ddc sample was taken         */
    double   ddc_Q;           /**< ddc Q value                                           */
    double   ddc_I;           /**< ddc I value                                           */
    double   ddc_amp;         /**< downconverted amplitude                               */
    double   ddc_phase;       /**< downconverted phase                                   */
    double   ddc_tdecay;      /**< downconverted decay time of waveform                  */

    double   ddc_pos;         /**< calculated position from ddc                          */
    double   ddc_slope;       /**< calculated slope from ddc                             */

    double   ddc_ct_amp;      /**< last measured calibration tone amplitude for this bpm */
    double   ddc_ct_phase;    /**< last measured calibration tone phase for this bpm     */

    // FIT parameters
    int      fit_success;     /**< do we have proper fit info ?                          */
    double   fit_Q;           /**< fit Q value                                           */
    double   fit_I;           /**< fit I value                                           */
    double   fit_amp;         /**< fitted amplitude                                      */
    double   fit_phase;       /**< fitted phase                                          */
    double   fit_freq;        /**< fitted frequency (MHz)                                */
    double   fit_tdecay;      /**< fitted decay time of waveform (usec)                  */
    double   fit_offset;      /**< fitted offset for waveform                            */

    double   fit_pos;         /**< calculated position from fit                          */
    double   fit_slope;       /**< calculated slope from fit                             */

    double   fit_ct_amp;      /**< last measured calibration tone amplitude for this bpm */
    double   fit_ct_phase;    /**< last measured calibration tone phase for this bpm     */
  };
  
  /**
     This structure contains the global beam parameters as well as 
     a pointer to the array of bunches
  */
  struct beamconf {
    int		  train_num;        /**< seq number of the train (evt num)               */

    double        beamrate;         /**< beam repetition rate (train to train)           */
    double        bunchrate;        /**< bunch repetition rate (in the train)            */
    int           nbunches;         /**< number of bunches per train                     */
    
    bunchconf_t  *bunch;            /**< list of pointers to the bunch conf structures   */

    double        position[2];      /**< beam position at the origin                     */
    double        positionsigma[2]; /**< position spread at the origin                   */

    double        slope[2];         /**< beam slope at the origin                        */
    double        slopesigma[2];    /**< slope spread at the origin                      */

    double        tilt[2];          /**< bunch tilt at the origin                        */
    double        tiltsigma[2];     /**< tilt spread at the origin                       */

    double        bunchlength;      /**< bunch length at the origin                      */
    double        bunchlengthsigma; /**< length spread at the origin                     */

    double        energy;           /**< beam energy (in GeV) at the origin              */
    double        energysigma;      /**< beam energy spread                              */
    double        charge;           /**< bunch charge (in nC)                            */
    double        chargesigma;      /**< charge spread                                   */
  };

  /**
     This structure contains information on a single bunch inside the 
     bunchtrain, which has its own energy, internal energy spread, charge,
     length, position/slope/tilt in the world coo frame and
     position/slope/tilt in the BPM local coo frame.
   */
  struct bunchconf {
    int      train_num;	      /**< seq number of the train this bunch belongs to         */
    int      bunch_num;	      /**< seq number of the bunch in the train                  */
    
    double   energy;          /**< energy of the bunch                                   */
    double   energyspread;    /**< energy spread inside the bunch                        */
    double   charge;
    double   length;          /**< the bunch length                                      */
    double   arrival_time;    /**< arrival time of bunch                                 */
    double   position[2];     /**< the bunch position x,y at the bpm coo                 */
    double   slope[2];        /**< the bunch slope x',y' at the bpm coo                  */
    double   tilt[2];         /**< the bunch tilt x',y' at the bpm coo                   */
    
    double   bpmposition[3];  /**< where the beam hits the BPM in the BPM local co       */
    double   bpmslope[2];     /**< slope of beam through the BPM in BPM local co         */   
    double   bpmtilt[2];      /**< bunch tilt in the BPM local co                        */
  };
  
  /**
     This structure defines a BPM resonant mode which is defined by 
     it's resonant frequency, Q factor and sensitivities to the
     beam charge, slope and bunch tilt. */
  struct bpmmode {
    char          name[20];     /**< The name for the BPM mode, e.g "dipolex" */
    double        frequency;    /**< The resonant frequency of the mode */
    double        Q;            /**< The Q factor for the mode */
    int           order;        /**< The mode order, 0:monopole, 1:dipole, 2:quadrupole... */
    enum bpmpol_t polarisation; /**< The mode polarisation: horiz, vert */
    double        sensitivity;  /**< The sensitivity of the mode, units depend on order */
    complexwf_t   *response;    /**< Pointer to the mode response buffer */
    doublewf_t    *buffer;      /**< Pointer to the mode's buffer */
  };

  /**
     This structure contains the complete RF model for a BPM, which is
     essentially a collection of it's resonant modes and sensitivities */
  struct rfmodel {
    char       name[20];        /**< A name for the cavity's RF model */
    int        nmodes;          /**< The number of BPM modes in the model */
    bpmmode_t *mode;            /**< A list of pointers to the array of modes */
  };




/* -----------------------------------------------------------------------------
   function prototypes and declarations
   -------------------------------------------------------------------------- */
  EXTERN int bpm_verbose;    /**< be a bit verbose in libbpm */
  EXTERN int libbpm_evtnum;  /**< the global event number in the processing */
  

#ifdef __cplusplus
}
#endif


#endif /* #ifndef BPMINTERFACE_H__ */
/*@}*/
/* ================================ end of file ============================= */
