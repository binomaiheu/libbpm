/**
   online_proc

   An example of an online processing algorithm that can be used at either KEK or
   ESA. The epics_post and epics_read routines are merely stubs to indicate where
   the proper calls would go.

   Author : Bino Maiheu, University College London (2006-2007)
            Mark Slater, University of Cambridge (2005-2007)

   ChangeLog :
*/

#include <pthread.h>

#include <bpm/bpm_units.h>
#include <bpm/bpm_interface.h>
#include <bpm/bpm_process.h>
#include <bpm/bpm_orbit.h>
#include <bpm/bpm_nr.h>
#include <bpm/bpm_calibration.h>
#include <bpm/bpm_analysis.h>

int gRunning;                   // Is the online processing running?
int gNumSIS;                    // The number of SISs we have
int gNumBPMs;                   // The number of BPMs
bpmsignal_t *gSignalCurrEvt;    // The current signal event
bpmproc_t *gProcCurrEvt;        // The current process event
bpmconf_t *gBPMConf;            // The BPM configurations
bpmcalib_t *gCalib;             // The BPM calibratons

// For resolution calcs
#define NUM_SVD 200             // Number of events to use in the SVD
bpmproc_t **gResBuf;            // The buffer for events
int gResBPMIdx[] = { 1, 0, 4 }; // The BPMs to use in the regression

// For calibrations
#define MAX_CAL 4000
#define CAL_FILENAME "examples/data/ESACalib.asc"
int gInCal = 0;                     // Are we in a calibration?
int gProcessingCal = 0;             // Are we processing a cal?
int gLoadCal = 0;                   // Should we load a cal?
int gCalEvt = 0;                    // The calibration event number
bpmsignal_t **gCalBuf;              // Buffer for calibrations
void *online_calibration (void *arg );  // The calibration routine


int main( int argc, char **argv ) {

  /* -----------------------------------------------------------
     The main entry point. Simply calls the other routines
     ------------------------------------------------------------
  */

  // Initialise routine
  init_online( );

  // Main loop
  main_loop( );

  // Clean up
  clean_up( );

}


int init_online( ) {

  /* -----------------------------------------------------------
     Initialise the various bits and pieces
     ------------------------------------------------------------
  */
  int ibpm;

  // We're running
  gRunning = 1;

  // Initialise EPICS
  epics_init( );

  // Read in initial stuff
  double sis;
  epics_read( "num_sis", &sis );
  gNumSIS = (int) sis;

  // Load in the bpm config info
  load_bpmconf( "examples/data/esa_bpm_conf.asc", &gBPMConf, &gNumBPMs );

  // Setup buffering
  gResBuf = (bpmproc_t **) calloc( 3, sizeof(bpmproc_t*) );   // 3 BPMs used for the residual
  for ( ibpm = 0; ibpm < 3; ibpm++) {
    gResBuf[ibpm] = (bpmproc_t *) calloc( NUM_SVD, sizeof(bpmproc_t) );
  }

  gCalBuf = (bpmsignal_t **) calloc( gNumBPMs, sizeof(bpmsignal_t*) );   // 3 BPMs used for the residual
  for ( ibpm = 0; ibpm < gNumBPMs; ibpm++) {
    gCalBuf[ibpm] = (bpmsignal_t *) calloc( MAX_CAL, sizeof(bpmsignal_t) );
  }

  // Set up event structures
  gSignalCurrEvt = (bpmsignal_t *) calloc( gNumBPMs, sizeof(bpmsignal_t*) );
  gProcCurrEvt = (bpmproc_t *) calloc( gNumBPMs, sizeof(bpmproc_t*) );

  for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
    gSignalCurrEvt[ibpm].wf = calloc( gBPMConf[ibpm].digi_nsamples, sizeof( int ) );
    gSignalCurrEvt[ibpm].ns = gBPMConf[ibpm].digi_nsamples;
  }

  // Load the calibration
  gCalib = (bpmcalib_t*) calloc( gNumBPMs, sizeof(bpmcalib_t) );
  load_calibration( CAL_FILENAME, gBPMConf, gCalib, gNumBPMs );

  // Post initial conditions to EPICS
  epics_post( "online_state", gRunning );
}

int main_loop( ) {

  /* -----------------------------------------------------------
     The main loop that keeps going as we recieve events
     ------------------------------------------------------------
  */

  int ibpm, ievt, in_cal = 0;
  char buf[255];

  while ( gRunning ) {

    //-----------------------------------------------------
    // Wait for trigger
    double trig = 0.0;
    while (!trig) {
      epics_read("trigger", &trig );
    }

    //-----------------------------------------------------
    // Check for the start or end of a calibration
    
    // This could be a magnetic check or taken from BPM info
    if (!gProcessingCal) {
     
      epics_read( "in_calibration", &in_cal );
      if ( ( in_cal ) && ( ! gInCal ) ) {
	
	// Initialise the calibration
	gCalEvt = 0;      
      }
      
      if ( ( !in_cal ) && ( gInCal ) ) {
	if (gCalEvt < MAX_CAL ) {
	  
	  // Fire off a thread to compute the calibraton
	  pthread_t tid;
	  void *arg;
	  pthread_create( &tid, NULL, online_calibration, arg );  
	  gProcessingCal = 1;
	} else {
	  
	  // The buffer's full and we haven't finished yet - abort!
	}
      }  
      
      gInCal = in_cal;
      if (gInCal) gCalEvt++;
    }

    // Should we load a new cal?
    if (gLoadCal)
      {
	load_calibration( CAL_FILENAME, gBPMConf, gCalib, gNumBPMs );
	gLoadCal = 0;
      }

    //-----------------------------------------------------
    // Grab the latest SIS data
    // (This assumes the bpm number is in the channel number order)
    for (ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {

	double val;
	int s;

	for (s = 0; s < gBPMConf[ibpm].digi_nsamples; s++ ) {

	  sprintf(buf, "sis_data[%d][%d][%d]", ibpm / 8, ibpm % 8, s);	
	  epics_read( buf, &val );
	  gSignalCurrEvt[ibpm].wf[s] = (int) val;
	}

    }

    // Store the waveforms in the calibration buffer if necessary
    if (gInCal && !gProcessingCal) {
      for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
	gCalBuf[ibpm][gCalEvt] = gSignalCurrEvt[ibpm];
	memcpy( gCalBuf[ibpm][gCalEvt].wf, gSignalCurrEvt[ibpm].wf, gCalBuf[ibpm][gCalEvt].ns * sizeof(int) );
      }
    }
      
    //-----------------------------------------------------
    // Now process the waveforms
    for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
      if ( gBPMConf[ibpm].cav_type == diode ) {
	process_diode( &(gBPMConf[ibpm]), &(gSignalCurrEvt[ibpm]), &(gProcCurrEvt[ibpm]) );	
      }
    }
    
    for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
      if ( gBPMConf[ibpm].cav_type == monopole ) {
	process_monopole(&(gBPMConf[ibpm]), &(gCalib[ibpm]),
			 &(gSignalCurrEvt[ibpm]), 
			 &(gProcCurrEvt[ibpm]),
			 &(gProcCurrEvt[ gBPMConf[ibpm].diode_idx ]),
			 PROC_DO_DDC);	
      }
    }
    
   for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
     if ( gBPMConf[ibpm].cav_type == dipole ) {
	process_dipole(&(gBPMConf[ibpm]), &(gCalib[ibpm]),
		       &(gSignalCurrEvt[ibpm]), 
		       &(gProcCurrEvt[ibpm]),
		       &(gProcCurrEvt[ gBPMConf[ibpm].diode_idx ]),
		       &(gProcCurrEvt[ gBPMConf[ibpm].ref_idx ]),
		       PROC_DO_DDC);
      }
    }   

   //-----------------------------------------------------
   // Now calculate resolutions, etc.
   
   // shift buffer and transfer the latest data
   for ( ibpm = 0; ibpm < 3; ibpm++ ) {
     for ( ievt = 0; ievt < NUM_SVD - 1; ievt++ ) {

       gResBuf[ibpm][ievt] = gResBuf[ibpm][ievt + 1];
     }

     gResBuf[ gResBPMIdx[ibpm] ][NUM_SVD - 1] = gProcCurrEvt[ibpm];
   }

   // Compute the SVD coefficients
   double coeffs[5];
   ana_get_svd_coeffs( gResBuf, 3, NUM_SVD, NUM_SVD, coeffs, ANA_SVD_TILT );

   // Find the resolution
   double mean, rms;
   ana_compute_residual( gResBuf, 3, NUM_SVD, coeffs, ANA_SVD_TILT, &mean, &rms );
   
   //-----------------------------------------------------
   // Finally, output all the processed info to EPICS
   for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
     
     sprintf(buf, "%s_pos", gBPMConf[ibpm].name );
     epics_post(buf, gProcCurrEvt[ibpm].ddc_pos );

     sprintf(buf, "%s_slope", gBPMConf[ibpm].name );
     epics_post(buf, gProcCurrEvt[ibpm].ddc_slope );

     epics_post("residual", rms );
   }

  } // while ( gRunning ) {
}
int clean_up( ) {

  /* -----------------------------------------------------------
     Free any allocated space, etc.
     ------------------------------------------------------------
  */
}

int epics_post( char* chan_name, double val ) {

  /* -----------------------------------------------------------
     Purely a stub
     ------------------------------------------------------------
  */

  return 1;
}

int epics_read( char* chan_name, double *val ) {

  /* -----------------------------------------------------------
     Purely a stub
     ------------------------------------------------------------
  */

  return 1;
}

int epics_init( char* chan_name, double val ) {

  /* -----------------------------------------------------------
     Purely a stub
     ------------------------------------------------------------
  */

  return 1;
}

void *online_calibration (void *arg )
{
  /* -----------------------------------------------------------
     Recalibrate the BPMs, store the result and report back
     ------------------------------------------------------------
  */

  int ievt, ibpm, num_events = gCalEvt;
  bpmproc_t cal_proc[gNumBPMs][num_events];
  bpmcalib_t cal_data[gNumBPMs];

  // Set up the initial calibrations
  for ( ibpm = 0; ibpm < gNumBPMs; ibpm++) {
    
    cal_data[ibpm].t0Offset = 0.2 * usec;
    cal_data[ibpm].ddcfiltBW = 4.* MHz;
    cal_data[ibpm].ddcepsFilt = 0.001;
  }

  //-----------------------------------------------------
  // First, loop over the events and fit them
  for ( ievt = 0; ievt < num_events; ievt++ )
    {

     // Process all bpms in this order: diode, reference, dipole
     for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
      if ( gBPMConf[ibpm].cav_type == diode ) {
	process_diode( &(gBPMConf[ibpm]), &(gCalBuf[ibpm][ievt]), &(cal_proc[ibpm][ievt]) );	
      }
    }
    
    for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
      if ( gBPMConf[ibpm].cav_type == monopole ) {
	process_monopole(&(gBPMConf[ibpm]), &(cal_data[ibpm]),
			 &(gCalBuf[ibpm][ievt]), 
			 &(cal_proc[ibpm][ievt]),
			 &(cal_proc[gBPMConf[ibpm].diode_idx][ievt]),
			 PROC_DO_FIT);	
      }
    }
    
   for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
     if ( gBPMConf[ibpm].cav_type == dipole ) {
	process_dipole(&(gBPMConf[ibpm]), &(cal_data[ibpm]),
		       &(gCalBuf[ibpm][ievt]),
		       &(cal_proc[ibpm][ievt]),
		       &(cal_proc[gBPMConf[ibpm].diode_idx][ievt]),
		       &(cal_proc[gBPMConf[ibpm].ref_idx][ievt]),
		       PROC_DO_FIT);
      }
    }   
   }

  //-----------------------------------------------------
  // Calculate the median frequencies and decay constants
  for (ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {

    if ( gBPMConf[ibpm].cav_type != diode ) {
      double freq[num_events], tdecay[num_events];
      
      for ( ievt = 0; ievt < num_events; ievt++ ) {

	// Transfer to temporary storage
	freq[ievt] = cal_proc[ibpm][ievt].fit_freq;
	tdecay[ievt] = cal_proc[ibpm][ievt].fit_tdecay;
      }
      
      cal_data[ibpm].freq = nr_median( num_events, freq );
      cal_data[ibpm].tdecay = nr_median( num_events, tdecay );
    }
  }
  
  //-----------------------------------------------------
  // Reprocess the waveforms to get positions
  for ( ievt = 0; ievt < num_events; ievt++ )
    {
      // refs
      for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
	if ( gBPMConf[ibpm].cav_type == monopole ) {
	  process_monopole(&(gBPMConf[ibpm]), &(cal_data[ibpm]),
			   &(gCalBuf[ibpm][ievt]), 
			   &(cal_proc[ibpm][ievt]),
			   &(cal_proc[gBPMConf[ibpm].diode_idx][ievt]),
			   PROC_DO_DDC);	
	}
      }
      
      for ( ibpm = 0; ibpm < gNumBPMs; ibpm++ ) {
	if ( gBPMConf[ibpm].cav_type == dipole ) {
	  process_dipole(&(gBPMConf[ibpm]), &(cal_data[ibpm]),
			 &(gCalBuf[ibpm][ievt]),
			 &(cal_proc[ibpm][ievt]),
			 &(cal_proc[gBPMConf[ibpm].diode_idx][ievt]),
			 &(cal_proc[gBPMConf[ibpm].ref_idx][ievt]),
			 PROC_DO_DDC);
	}
      }   
    }
  
  //-----------------------------------------------------
  // Finally, calibrate the BPMs
  beamconf_t calc_beampos[gNumBPMs][num_events];
  int cal_steps = 5;
 
  // *****************************************************************************
  // This is where the expected bema position is calculated. Change this for the
  // approrpiate corrector/mover scan values
  // *****************************************************************************

  for (ibpm = 0; ibpm < gNumBPMs; ibpm++)
    {
      if ( ( gBPMConf[ibpm].cav_type == dipole ) && ( gBPMConf[ibpm].cav_polarisation == horiz ) ) {

	setup_calibration( &(gBPMConf[ibpm]), &(cal_proc[ibpm][0]), num_events,
			   0, num_events, 0. * urad, -170., 
			   170., cal_steps, &(calc_beampos[ibpm][0]) );
	
	calibrate( &(gBPMConf[ibpm]), &(calc_beampos[ibpm][0]), &(cal_proc[ibpm][0]), 
		   num_events, &(cal_data[ibpm]) );
      }
    }
  
  save_calibration( CAL_FILENAME, gBPMConf, cal_data, gNumBPMs );

  // We've finished, so get the main thread to load up the new cal when ready
  gProcessingCal = 0;
  gLoadCal = 1;

}
