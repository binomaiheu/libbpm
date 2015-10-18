/*
  libbpm - BPM signal processing/simulation library
  Copyright (C) 2006-08 Bino Maiheu (bino@hep.ucl.ac.uk)
  
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
    @defgroup processing BPM Processing Routines
    @anchor proc 

    This set of routines contains the BPM digitised waveform processing routines
    to go from a sis digitised waveform to position and slope information. 

    @section proc_structure General structure of the BPM signal processing

    The BPM signal processing algorithms are centered around a few top-level routines which 
    need to called by a standard user. All make use of a number of BPM data structures
    which hold BPM configuration data ( bpmconf_t ), processed BPM information ( bpmproc_t )
    or BPM calibration information ( bpmcalib_t ). As the BPM processing algorithms make
    extensive use of the bpmdsp module, the BPM signals need to be encapsulated in a 
    doublewf_t waveform before feeding them to these processing routines. The top-level
    processing routines have a mode bitword which provides some processing options that
    the user can feed into the processing algorithm. 

    @subsection proc_structure_diode Diode signal processing

    Since the idea was to unify the processing into one coherent set of data structures, 
    the diode or trigger information had to be fitted into the same framework as the
    BPM data. This is the function call : 

    @code
    int process_diode( doublewf_t *signal, bpmconf_t *conf, bpmproc_t *proc );
    @endcode

    So the diode pulse has to be fitted into a doublewf_t along with a bpmconf_t structure
    conf. The routine first checks the flag bpmconf_t::cav_type for the cavity type. This
    should be of type diode for the routine to proceed. It then calls the fit_diodepulse
    routine onto the signal, which returns the fitted t0 into the bpmproc_t structure as
    proc->t0. 

    @attention
    Note that there is the possibility to abuse a dipole or monopole signal as a trigger
    pulse. In this case the process_diode routine will determine the RMS of the noise 
    in front of the digitised dipole/monopole signal (first 20 samples ) and return the
    timestamp in bpmproc_t::t0 of the first sample which is 10 times largers than this 
    RMS value. For this behaviour, the bpmconf_t::cav_type setting is irrelevant but the
    bpmconf_t::forced_trigger value has to be set to 1. Note that this behaviour is normally
    not needed an for experimental purposes only.
    @endattention


    @subsection proc_structure_monopole Monopole signal processing
    
    For monopole cavities one only needs to determine the amplitude and phase, so 
    no post-processing to get to position and slope using a reference cavity and
    calibration information is needed. Therefore the process_monopole routine
    is basically a wrapper around the process_waveform routine which does
    exactly this determination of the amplitude and phase. The function call is :

    @code
    int process_monopole( doublewf_t *signal, bpmconf_t *bpm, bpmproc_t *proc, 
                          bpmproc_t *trig, unsigned int mode );
    @endcode
    
    This routine basically is a wrapper around

    @code
    int process_waveform( doublewf_t *signal, bpmconf_t *bpm, bpmproc_t *proc, 
                          bpmproc_t *trig, unsigned int mode );
    @endcode

    and handles all the processing steps flagged by the mode bitword. Chronologically
    it executes the following steps :

    -  Check whether the waveform was saturated or not. This is done by a call to 
       check_saturation, which needs the doublewf_t signal obviously and the ADC
       resolution set by the number of bits in bpmconf_t::digi_nbits. It returns 
       whether the waveform was saturated ( saved in bpmproc_t::saturated ) and
       assigns the sample number of the first unsaturated sample in the waveform
       to bpmproc_t::iunsat.

    -  Then process_waveform goes on with subtracting the pedestal of the waveform
       by getting the average and RMS of the first 20 samples in the waveform using
       get_pedestal and storing the results in bpmproc_t::voltageoffset and 
       bpmproc_t::ampnoise. It subsequently subtracts this voltage offset from each
       sample in the waveform. 
    
    -  Then the t0 time is set. If the process_waveform has trigger information
       available in the form of a bpmproc_t trigger argument wich was handled by 
       process_diode, then the routine will assume this information has to be used
       as t0 and will copy the trigger->t0 value to it's own bpmproc_t::t0. If a 
       bpmproc_t trigger argument is not available ( NULL pointer ), the 
       process_waveform routine will assume the t0 has been set fixed by the BPM
       configuration ( external clocking ) and will use and copy the bpmconf_t::t0
       to it's bpmproc_t::t0 value. The bpmproc_t::t0 is furtheron used in the rest
       of the processing as the starting time for this cavity signal. 

    -  If the PROC_DO_FFT flag has been set in the mode bitword, the process_waveform
       routine will compute the waveform FFT by calling fft_waveform from the bpmdsp
       module and storing the result in bpmproc_t::ft. If this is succesfull, the 
       code will go on to check whether this fourier transform needs to be fitted for 
       it's frequency and decay time ( Lorentz line width ). This is done by calling
       fit_fft. 
       @attention
       This routine is a little experimental and can easily by replaced by the user
       with some other package e.g. ROOT. The full complex fourier waveform is available
       in the bpmproc_t::ft as a complexwf_t.
       @endattention

    -  If the PROC_DO_FIT flag has been set in the mode bitword, the process_waveform
       routine will try to fit a decaying sinewave to the waveform, attempting to extract
       amplitude, phase, frequency and decay time. 
       \attention This routine is quite experimental as well and needs proper checking 
                  before it can be used stabily ! I recommend using a proper fitting 
		  package such as MINUIT to fit the waveforms to a decaying sine wave. 


    -  If the PROC_DO_DDC flag has been set in the mode bitword, the process_waveform
       routine will perform the digital downconversion on the waveform. As this is 
       a more complex algorithm, we will go into a bit more detail here. 
       - First, we have to tell the DDC algoritm where to get it's frequency and 
         decay time from. By default the algorithm will use in both cases the
         frequency and decay time which are set in the cavities configuration, being
         bpmconf_t::ddc_freq and bpmconf_t::ddc_tdecay. However, if the flag(s)
         PROC_DDC_FITFREQ and/or PROC_DDC_FITTDECAY is/are present and the
	 fits ( see previous item ) were succesfull, the ddc algorithm will use the 
	 fitted frequency and decaytime values. Alternatively, if the flag(s) 
         PROC_DDC_FFTFREQ and/or PROC_DDC_FFTTDECAY are/is present, the ddc algoritm 
         will use the frequency and decay time derived from the fitted lorentz lineshape 
         of the waveforms fourier transform. 

       - Next the DDC algoritm handls the saturation if present ( was set by the 
         bpmproc_t::saturated ) flag already. If the waveform was saturated, we will
	 shift the position of the sample time to the last unsaturated sample. 
	 \attention Since people haven't converged on a proper way to handle saturation, 
                    this is a bit of an open point in the code. At the moment, the 
                    ddc_tSample is set to the last unsaturated sample, but one should 
                    take into account somehow the bandwidth of the DDC filter, which is 
                    not done. I've left it as it is, with the wise advice to store the 
                    bpmproc_t::saturated flag into the user data and simply cut away 
                    those pulses. 

	 If no saturation is present, the sampling point (expressed in time-units, 
	 not sampled ) of the DDC algoritm is set to the t0 time ( starting point of 
	 the waveform ) + a constant time offset, which can be tweaked in optimisation.
	 @code
	 proc->ddc_tSample = proc->t0 + bpm->ddc_tOffset;
	 @endcode

       - After the sampling time has been calculated in the previous step, it is converted
         into a sample number and stored in bpmproc_t::ddc_iSample. 

       - Then the real downconversion is done, by default libbpm will try to use the
         optimised ddc_sample_waveform routine to save CPU cycles, but if the
         full DDC is requested by the mode flag PROC_DDC_FULL, it will go through the 
	 entire waveform and convert it to DC using the frequency set as explained previously. 
         The routine that is called is ddc_waveform which basically needs the the 
	 pedestal subtraced doublewf_t waveform, the frequency of downconversion, a 
	 2 omega filter, defined by a filter_t structure having the correct type (lowpass) 
         and bandwidth already define and stored in bpmconf_t::ddc_filter. The full complex 
         downconverted waveform is stored in the case of full ddc in bpmproc_t::dc. The 
         amplitude and phase are calculated at the t0 time by extrapolating the phase and 
	 amplitude back from the sampling point at bpmproc_t::ddc_iSample. The 
	 ddc_sample_waveform returns these values directly, but does it internally by 
	 extrapolation from the sampling time as well, one therefore needs to provide t0, 
	 tdecay and iSample as additional arguments to ddc_sample_waveform compared to 
	 ddc_waveform. 

       - After this is done, the determined phase is normalised in between 0 and 2pi. 

    @subsection proc_structure_dipole Dipole signal processing
    
    Dipole cavity waveforms first need to undergo the same processing step as monopole 
    waveforms, to determine their phase and amplitude. After that position and 
    slope information need to be determined using the calibration information. The routine
    
    @code
    int process_dipole( doublewf_t *signal, bpmconf_t *bpm, bpmcalib_t *cal, bpmproc_t *proc, 
     	                bpmproc_t *trig, bpmproc_t *ampref, bpmproc_t *phaseref, 
			unsigned int mode );
    @endcode

    is therefore a wrapper around the following two core routines :

    @code
    int process_waveform( signal, bpm, proc, trig, mode );
    int postprocess_waveform( bpm, proc, cal, ampref, phaseref, mode );
    @endcode

    \attention If the PROC_CORR_GAIN (or PROC_CORR_AMP, PROC_CORR_PHASE) flag is set in the 
    mode word, the process_dipole routine will correct the gains based upon the latest 
    calibration tone information stored in the bpmproc_t::ddc_ct_amp etc variables and 
    comparing them to the bpmcalib_t::ddc_ct_amp at the time of calibration. This is done 
    by a call to correct_gain

    The process_waveform is explained under the process_monopole cavity, the 
    postprocess_waveform routine executes the following :

    - Firstly the routine calculates the I and Q for the dipole cavity from the
      amplitude and phase references. This is done by a call to get_IQ, and the
      values are stored in bpmproc_t::ddc_Q, bpmproc_t::ddc_I for the DDC information
      and bpmproc_t::fit_Q, bpmproc_t::fit_I for the fitted information.

    - For dipole cavities, the real phase information that means anything is the
      phase difference between the reference cavity and the dipole cavity. This get's 
      stored into bpmproc_t::ddc_phase and/or bpmproc_t::fit_phase. If the flag
      PROC_RAW_PHASE is set in the mode word, this is skipped. 

    - Using the I and Q information, the position and slope are calculated 


    @section proc_flow Processing flow
    The question now is how to organise the processing flow from the digitised waveform
    data. Before being able to obtain positions and slopes, the user will need to have
    processed all the trigger ( diode ) pulses. And thereafter the monopole waveforms 
    in the event. After that positions and slopes can be calculted using the process_diopole 
    routine. Note that the monopole waveforms depend on the trigger information in the
    case of internal triggering using a trigger pulse, so a good way to proceed is 
    first to all the trigger pulses, than all the monopole pulses and then all the dipole
    waveforms.

    Alternatively the user can first use the routine process_waveform on all of the 
    waveforms ( together with processing the trigger information ). After this is
    done, the user can use the postprocess_waveform routine to perform the post-processing
    on the dipole waveforms. 

    @section proc_timing About trigger pulses, internal vs. external clock
    
    The SIS ADCs can be triggered by using an external clock in which case all the
    modules in the system are synchronised and no trigger pulses are needed. Because of the
    way the processing is setup in process_waveform, the user has to be mindfull of a 
    number of things depending on whether the ADC modules are triggered internally ( and
    a trigger pulse is available ) or whether they are triggered externally, synchronised
    to the beam clock, in which case the starting time ( t0 ) of the pulses should be
    constant for each individual BPM signal. 

    @subsection proc_timing_ext External clock triggering

    In this case, the t0 should be set in the BPM configuration under bpmconf_t::t0. During
    the processing this value will be used and copied to bpmproc_t::t0. The bpmconf_t::tOffset
    defines the offset from this t0 of the pulse of the sampling point in the waveform such 
    that 
    @code
    proc->ddc_tSample = proc->t0 + bpm->ddc_tOffset;
    @endcode
    This mode will be assumed automatically in the absence of the 4th argument of 
    process_waveform ( bpmproc_t *trig = NULL ). 

    @subsection proc_timing_int Internal clock triggering
    
    There the bpmconf_t::t0 value is ignored and no t0 value needs to be specified before
    hand since it will be fitted from the diode/trigger pulse. In this case the 4th argument
    of process_waveform needs to be present. Also,<b>the bpmconf_t::tOffset keeps it's
    definition exaclty the same as in the external clock case</b>. It is the time difference
    between the sample time and the starttime of the waveform t0, which in this case
    got fit instead of being fixed.


    @section Incorporating calibration tone information

    The calibration tone information is kept in two locations. Firstly at the time of
    calibration, the user should make sure that the latest calibration tone information
    is set in the bpmcalib_t structure under bpmcalib_t::ddc_ct_amp and bpmcalib_t::ddc_ct_phase
    and analoguous for the parameters for the fitted processing. Than each time a 
    calibration tone pulse is encountered, the user should pass the phase and amplitude 
    of the calibration tone on to the bpmproc_t::ddc_ct_amp and bpmproc_t::ddc_ct_phase and 
    therefore always keep the lateste calibration tone information in this location. 
    Each call to 
    
    @code
    int correct_gain( bpmproc_t *proc, bpmcalib_t *cal, unsigned int mode )
    @endcode

    then corrects the phase and amplitude of the current pulse by scaling the amplitude
    with the ratio between the caltone amplitude at the time of calibration and the 
    lastest one and shifting the phase by the phase difference between the phase of the
    calibration tone at the time of BPM calibration and the latest phase recorded in
    the bpmproc_t::ddc_ct_phase variable ( or bpmproc_t::fit_ct_phase ).

    \attention I've include a mode bitword, which takes the flags PROC_CORR_GAIN to correct
    both amplitude and phase, and PROC_CORR_AMP, PROC_CORR_PHASE to correct only 
    one parameter individually. This is done since e.g. for internal clocking, when
    the ADC's are not synchronised to each other, it is not really clear where to sample
    the waveform unless a trigger is supplied in the ADC. For external synchronized
    clocking, we can just give a fixed sample number, stored in the bpm configuration
    under bpmconf_t::ddc_ct_iSample.

*/    


/** 
    @file  
    @ingroup processing
    @brief   libbpm main processing routines
    
    This header contains the definitions for libbpm's main BPM processing routines
*/

/** @addtogroup processing */
/** @{ */

#ifndef BPMPROCESS_H__
#define BPMPROCESS_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <float.h>
#include <math.h>
#include <bpm/bpm_defs.h>
#include <bpm/bpm_interface.h>
#include <bpm/bpm_wf.h>
#include <bpm/bpm_dsp.h>

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */
#define  PROC_DEFAULT         0           //! default mode
#define  PROC_DO_FFT          0x00000001  //! do fft + fit for fft_freq and fft_tdecay
#define  PROC_DO_FIT          0x00000002  //! do fit, gets automatically freq + tdecay
#define  PROC_DO_DDC          0x00000004  //! do ddc, needs to know where to get freq + tdecay
#define  PROC_DDC_CALIBFREQ   0x00000008  //! ddc gets freq from calib structure
#define  PROC_DDC_CALIBTDECAY 0x00000010  //! ddc gets freq from calib structure
#define  PROC_DDC_FITFREQ     0x00000020  //! ddc gets freq from fit
#define  PROC_DDC_FITTDECAY   0x00000040  //! ddc gets tdecay from fit
#define  PROC_DDC_FFTFREQ     0x00000080  //! ddc gets freq from fft
#define  PROC_DDC_FFTTDECAY   0x00000100  //! ddc gets freq from tdecay
#define  PROC_DDC_FULL        0x00000200  //! do the full DDC of the waveform
#define  PROC_FIT_DDC         0x00000400  //! fit the ddc for tdecay
#define  PROC_FIT_FFT         0x00000800  //! fit the fft using libbpm
#define  PROC_RAW_PHASE       0x00001000  //! leave the phase uncorrected in bpmproc_t
#define  PROC_CORR_AMP        0x00002000  //! correct the amplitude using caltone information
#define  PROC_CORR_PHASE      0x00004000  //! correct the phase using caltone information
#define  PROC_CORR_GAIN       0x00006000  //! correct both amplitude and phase

/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif


/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */
  /* -- General processing routines to process monopole, dipole and diode pulses -- */

  /**
     This routine processes a diode pulse, which should be found in the
     signal structure. It fills the proc structure with the t0. The routine 
     checks what the signal type (conf->cav_type) is and when it really is a 
     diode pulse, it will fit the pulse and return t0, otherwise (when the 
     signal is a monopole or dipole signal), it will determine the onset 
     of the waveform by looking where the signal's absolute value exceeds 
     10 * the noise RMS at the beginning of the waveform.
     
     @param signal  The bpm signal
     @param conf  The bpm configuration structure
     @param proc The processed trigger structure (containing the t0)
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int process_diode( doublewf_t *signal, bpmconf_t *conf, bpmproc_t *proc );


  /**
     Top-level routine which is basically a wrapper around process_waveform and
     correct_gain to take into account the calibration tone data. See more in details
     documentation in those routines.
     
     @param signal The doublewf_t encoded BPM signal
     @param bpm    The bpm configuration structure
     @param cal    The bpm calibration structure, needed for the gain correction
     @param proc   The processed data structure
     @param trig   The structure with processed trigger info for that waveform
     @param mode   A bitpattern encoding what exactly to process
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int process_monopole( doublewf_t *signal, bpmconf_t *bpm, bpmcalib_t *cal, bpmproc_t *proc,
			       bpmproc_t *trig, unsigned int mode );

  
  /**
     Top-level routine which is a wrapper around process_waveform, correct_gain and
     postprocess_waveform. See more details in the documentation of those individual
     routines.
   
     @param signal   The doublewf_t encoded BPM signal
     @param bpm      The bpm configuration structure
     @param cal      The bpm calibration structure, needed for the gain correction
     @param proc     The processed data structure
     @param trig     The structure with processed trigger info for that waveform
     @param ampref   The already processed amplitude reference bpmproc_t structure
     @param phaseref The already processed phase reference bpmproc_t structure
     @param mode     A bitpattern encoding what exactly to process
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int process_dipole( doublewf_t *signal, bpmconf_t *bpm, bpmcalib_t *cal, bpmproc_t *proc, 
			     bpmproc_t *trig, bpmproc_t *ampref, bpmproc_t *phaseref, 
			     unsigned int mode );

  /**
     Top-level routine to processes a BPM beam pulse waveform (decaying "sin"-like wave) and
     derive amplitude and phase from the signal. The routine needs to be fed with a doublewf_t
     containing the digitized signal. The signal is checked for saturation, it's pedestal is
     determined and removed, the pulse starttime (t0) is set from the configuration or the
     trigger. Then, depending on the mode bitpattern, an FFT is performed, the waveform 
     is fitted and a digital downconversion is done. The results (amplitude and phase) are 
     stored in the bpmproc_t structure of the BPM.
     
     Relevant mode bit patterns for this routine are : 
     - PROC_DO_FFT : The Fourier Transform of the waveform gets computed and stored as a 
                     complexwf_t in the bpmproc_t::ft variable. 
     - PROC_FIT_FFT : An attempt to fit the Fourier Transform is made using a Lorentizan Lineshape.
                      If successfull, the bpmproc_t::fft_freq and bpmproc_t::fft_tdecay variables
		      will contain the fitted frequency and decaytime. I recommend however to
                      use a 3th party fitting routine for this (e.g. MINUIT) and implement this
                      in a user program.
     - PROC_DO_FIT : Attempts to fit a decaying sine wave to the waveform having the frequency,
                     the decay time, the amplitude and phase as free parameters. If successfull,
	  	     the bpmproc_t::fit_freq, bpmproc_t::fit_amp, bpmproc_t::fit_phase and 
		     bpmproc_t::fit_tdecay will contain the fit parameters. Again, I recommend
		     to use a 3th party fitting routine for this. 
     - PROC_DO_DDC : Will perform a digital downconversion on the waveform. The results are 
                     contained in bpmproc_t::ddc_amp and bpmproc_t::ddc_phase, determined at 
		     bpmproc_t::ddc_tSample, but extrapolated back to bpmproc_t::t0.
     - PROC_DDC_FITTDECAY, PROC_DDC_FFTTDECAY : Normally the ddc algoritm gets it's decay time for
                       extrapolation back to t0 from the bpmconf_t::ddc_tdecay variable, if one of
	  	       these flags are set it will get them from the fitted waveform or FFT if
                       they were succesful.
     - PROC_DDC_FITFREQ, PROC_DDC_FFTFREQ : Analogous as the previous item, but now for the ddc
                       frequency which is normally obtained from bpmconf_t::ddc_freq.
     - PROC_DDC_FULL : Will perform the DDC algorithm on the entire waveform and store the result
                       in bpmproc_t::dc
   
     @param signal The digitized signal converted into a doublewf_t
     @param bpm    A pointer to the bpmconf_t structure for the BPM channel
     @param proc   A pointer to the bpmproc_t structure for the BPM channel
     @param trig   A pointer to the bpmproc_t structure of the trigger for this BPM channel,
                   if this parameter is NULL, externall clocking will be assumed and the t0
                   from the bpmconf_t structure will be used in the processing. 
     @param mode   The processing mode bitword
     @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure 
  */
  EXTERN int process_waveform( doublewf_t *signal, bpmconf_t *bpm, bpmproc_t *proc, 
			       bpmproc_t *trig, unsigned int mode );
  



  /**
     Top-level routine to Post-process a waveform for whith the amplitude and the phase
     have already been defined using process_waveform. This routine goes on to calculate 
     I and Q from the phase and amplitudes as well as the postion and slope using the
     calibration information. 

     Relevant mode bit patterns for this routine are : 
     - PROC_RAW_PHASE : when this bit is active in the mode word, the routine will not 
                        replace the phase in the bpmproc_t structure by the phase difference
                        between the reference cavity and the processed cavity. Under normal
                        circumstances you don't want this since it's only the phase
                        difference which actually has any physical meaning.

     @param signal The digitized signal converted into a doublewf_t
     @param bpm    A pointer to the bpmconf_t structure for the BPM channel
     @param proc   A pointer to the bpmproc_t structure for the BPM channel
     @param cal    A pointer to the bpmcalib_t structure for the BPM channel
     @param ampref   A pointer to the bpmproc_t structure of the amplitude reference channel
                     for this BPM.
     @param phaseref A pointer to the bpmproc_t structure of the phase reference channel
                     for this BPM.
     @param mode   The processing mode bitword
     @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure
  */ 
  EXTERN int postprocess_waveform(  bpmconf_t *bpm, bpmproc_t *proc, bpmcalib_t *cal, 
				    bpmproc_t *ampref, bpmproc_t *phaseref, unsigned int mode );


  
  /* ===========================  caltone handling ===================================== */


  /**
     Top level routine to process the calibration tone via DDC, similar to process_waveform
     but it also updates the ddc_ct_amp and ddc_ct_phase variables in the bpmproc_t structure.
     No fitting is implemented in this routine.

     Relevant mode bit patterns for this routine are analogous as in process_waveform
     - PROC_DO_FFT : see process_waveform
     - PROC_FIT_FFT : see process_waveform
     - PROC_DO_DDC : see process_waveform

     @param signal The digitized signal converted into a doublewf_t
     @param bpm    A pointer to the bpmconf_t structure for the BPM channel
     @param proc   A pointer to the bpmproc_t structure for the BPM channel
     @param mode   The processing mode bitword

     @return BPM_SUCCESS upon succes, BPM_FAILURE upon failure
  */
  EXTERN int process_caltone( doublewf_t *signal, bpmconf_t *bpm, bpmproc_t *proc, 
			      unsigned int mode );


  /**
     Correct the processed amplitude and phase by using calibration tone information
     if the ddc and or fits were successfull. Since e.g. for internal clock it is 
     not really sure the phase information can be used if there is no proper trigger,
     some mode bits can be flagged to only correct the amplitude. 

     Relevant mode bit patterns for this routine are :
     - PROC_CORR_AMP : Correct the amplitude
     - PROC_CORR_PHASE : Correct the phase
     - PROC_CORR_GAIN : Correct both of them
     
     @param proc The bpmproc_t structure of the bpm
     @param cal  The bpmcalib_t structure of the bpm
     @param mode Mode of correction

     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int correct_gain( bpmproc_t *proc, bpmcalib_t *cal, unsigned int mode );



  /* ============================= Signal fitting routines ================================= */
  /**
     Fits the waveform with a decaying sin wave using the lmder/lmdif routines
     from nr_levmar.c !

     \attention Note that this routine is highly experimental, so don't use it for real
     production stuff. Instead I recommend using a proper minimisation package like
     MINUIT or so...
     
     @param *w      The waveform encoded as a doublewf_t
     @param  t0     t0 for the waveform
     @param i_freq    Initial frequency for the fit
     @param i_tdecay  Initial decay time for the fit
     @param i_amp     Initial amplitude for the fit
     @param i_phase   Initial phase for the fit
     @param freq      Fitted frequency
     @param tdecay    Fitted decay time
     @param amp       Fitted amplitude 
     @param phase     Fitted phase
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int fit_waveform( doublewf_t *w, double t0,
			   double i_freq, double i_tdecay, double i_amp, double i_phase,
			   double *freq, double *tdecay, double *amp, double *phase );

  
  /**
     Fits the diode pulse, basically a wrapper for get_t0, to conserve names and consistency
     in the library... is nothing more than a wrapper around get_t0, so see there...
  */
  EXTERN int fit_diodepulse( doublewf_t *w, double *t0 );


  /* ============================== FFT routines ====================================== */
  /**
     Performs a fast fourier transform of the waveform, after subtracting the 
     pedestal, basically just a wrapper around the forward realfft routine from 
     the DSP module. Please see it's documentation for more details...
  
     @param *w   the waveform
     @param  fft the complex returned fft spectrum
  
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int fft_waveform( doublewf_t *w, complexwf_t *ft );

  
  /**
     This routine prepares the fft fit of the waveform. It starts by getting the 
     position of the maximum in the spectrum (first nyquist band only). Then from
     this position runs left and right to determine where the amplitude drops to half
     of the peak amplitude and have an initial estimation of the FWHM. It will then
     set twice the FWHM width as the fit range in which to perform the fit, this is
     than returned by the samplnumbers n1 and n2.

     @param ft   The complexwf_t fourier transform 
     @param n1   The first sample to start the fit from
     @param n2   The last sample to take into account in the following fit
     @param amp  Initial estimation of the amplitude for the fit
     @param freq Initial estimation of the frequency for the fit
     @param fwhm Initial estimation of the FWHM for the fit.

     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure.
  */
  EXTERN int fit_fft_prepare( complexwf_t *ft, int *n1, int *n2,
			      double *amp, double *freq, double *fwhm );

  
  /**
     Fits the power spectrum of the FT of a waveform frequency and decay time.
     Internally it makes a call to fit_fft_prepare to get an initial estimation
     of the parameters and goes on by applying the nr_lmder routine to minimise
     the fourier transform power spectrum agains a lorentzian lineshape defined by 

     @f[ L = \frac{p_0}{ ( f - p_1 )^2 + \left( \frac{p_2}{2} \right)^2 } + p_3  @f]

     Where
     - p0 = the amplitude of the power spectrum
     - p1 = the frequency of the fourier transform peak
     - p2 = the full width at half maximum
     - p3 = a constant offset 
  
     @param  ft      The complexwf_t encoded fourier transform
     @param  freq    The returned frequency (p1)
     @param  tdecay  The returned tdecay (p2)
     @param  a       p0 (amplitude of powerspectrum ) of the fit ( can be NULL if not interested )
     @param  c       p3 (offset) of the fit ( can be NULL if not interested )
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int fit_fft( complexwf_t *ft, double *freq, double *tdecay, double *A, double *C );



  /* ============================== DDC routines ====================================== */

  /**
     Checks the saturation, so computes the first sample where no saturation
     occurs. If no saturation occurred in the waveform, this sample - stored in 
     iunsat - will be set to 0. 
     A saturated   sample is found when it's ADC value is more (resp. less) than then maximum
     allowed ADC value ( 2^nbits ) minus a threshold set to 15. ( resp. the minium
     allowed ADC value, being 0 ) plus a threshold set to 15. 
     
     \attention The waveform contained in the doublewf_t SHOULD NOT have been pedestal
     corrected. This routine will assume the waveform runs between 0 and 2^nbits. 

     Note the return code of the routine is slightly different than whan is conventional
     in libbpm since I wanted to encode whether saturation was found or not as the
     return code of the routine.

     @param w      The waveform to check, encoded as a doublewf_t
     @param nbits  The number of digitiser bits (e.g. 12 or 14 )
     @param iunsat The returned last unsaturated sample
   
     @return 1 when saturation was present, 0 when not, -1 when failure occurred
  */
  EXTERN int check_saturation( doublewf_t *w, int nbits, int *iunsat );


  /**
     Downmixes the input waveform agains a complex LO using a frequency f and phase 0, 
     the real part of the resulting complex waveform was mixed against a cosine-like
     wave, the imaginary part against a sinus-like. Note that this is just the downmixing
     itself, no filtering whatsoever is applied here. 
     
     @param w    The input waveform, encoded as a doublewf_t
     @param freq The frequency of the digital LO
     @param out  The complex output downmixed waveform
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure.
  */
  EXTERN int downmix_waveform( doublewf_t *w, double frequency, complexwf_t *out );


  /**
     As this is a pure wrapper around the ddc routine out of the dsp packate, please
     see the documentation there.
  */
  EXTERN int ddc_waveform( doublewf_t *w, double frequency, filter_t *filt, complexwf_t *dc,
			   doublewf_t *buf_re, doublewf_t *buf_im );

  
  /**
     TO BE IMPLEMENTED !!!

     This routine will contain a quicker version of the ddc algorithm that doesn't filter
     the entire waveform and only applies the filter at the sampling point. However, 
     I need to make custom a apply_filter routine which is universally valid for all types
     of filters (IIR as well).
  */
  EXTERN int ddc_sample_waveform( doublewf_t *w, double frequency, filter_t *filt, int iSample,
				  double t0, double tdecay, double *amp, double *phase,
				  doublewf_t *buf_re, doublewf_t *buf_im );

  /**
     Find the mean pedestal using the first 20 (or how ever many are required)
     sample values, store the results in the offset and rms. This routine in fact just
     calls the doublewf_basic_stats routine and gets the appropriate values from the
     wfstat_t structure.
     
     @param wf      The signal encoded as a doublewf_t
     @param range   The maximum sample to go to average over. The pedestal gets determined 
                    from the first "range" samples of the waveform
     @param *offset Returns the mean value of the samples, so voltage offset (pedestal value)
     @param *rms    Returns the RMS on that
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int get_pedestal( doublewf_t *wf, int range, double *offset, double *rms );

 
  /**
     Finds the t0 value from a diode peak, used in the case of internall triggering
     when a trigger pulse needs to be specified to calculate beam arrival

     \attention This routine needs some optimisation in terms of speed and some
     general checking in terms of correctness. Probably some re-writing using 
     the bpmwf structures would be good...
     
     @param  w   A pointer to the doublewf_t signal
     @param  t0  returns t0
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */ 
  EXTERN int get_t0( doublewf_t *w, double *t0 );
  

  
  /* ======================= Some position, IQ routines... ============================ */

  /**
     Gets the I and Q from the amplitude and phase of the waveform and it's respective
     references. The I and Q are calculated respectively as :
     @f[ I = \frac{A}{A_{ref}} \cos( \phi - \phi_{ref} )  @f]
     and
     @f[ Q = \frac{A}{A_{ref}} \sin( \phi - \phi_{ref} )  @f]
     
     @param amp      The amplitude of the considered waveform
     @param phase    The phase of the considered waveform
     @param refamp   The amplitude of the reference cavity
     @param refphase The phase of the reference cavity
     @param Q        The returned Q value
     @param I        The returned I value
     
     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int get_IQ( double amp, double phase, double refamp, double refphase,
		     double *Q, double *I );


  /**
     Returns the beam given I and Q values, IQphase and scale, it is calcualted as

     @f[ x = c \left[ I \cos (\phi_{IQ} ) + Q \sin ( \phi_IQ )\right] @f]
     
     Where c is the positionscale and x the position. 

     @param Q        The Q value (obtained from get_IQ)
     @param I        The I value (obtained from get_IQ)
     @param IQphase  The IQ phase rotation
     @param posscale The position scale 
     @param pos      The returned position

     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int get_pos( double Q, double I, double IQphase, double posscale,
		      double *pos );


  
  /**
     Returns the beam slope given I and Q values, IQphase and scale, it is calcualted as

     @f[ x' = c \left[ - I \sin (\phi_{IQ} ) + Q \cos ( \phi_IQ )\right] @f]
     
     Where c is the positionscale and x the position. 

     @param Q        The Q value (obtained from get_IQ)
     @param I        The I value (obtained from get_IQ)
     @param IQphase  The IQ phase rotation
     @param slopescale The slope scale 
     @param slope      The returned slope

     @return BPM_SUCCESS upon success, BPM_FAILURE upon failure
  */
  EXTERN int get_slope( double Q, double I, double IQphase, double slopescale,
			double *slope );

#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMINTERFACE_H__ */
/*@}*/
/* ================================ end of file ============================= */
