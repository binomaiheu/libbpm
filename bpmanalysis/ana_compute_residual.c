/**
   @file
   @ingroup analysis
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_analysis.h>

int ana_compute_residual( bpmproc_t **proc, int num_bpms, int num_evts,
			  double *coeffs, int mode, double *mean, double *rms ) {


  double *good_res = (double*) calloc( num_evts, sizeof(double));
  int evt, ibpm, idx, count = 0, good_evt;
  double residual;

  for ( evt = 0; evt < num_evts; evt++ ) {
    
    idx = 0;
    good_evt = 1;
    residual = proc[0][evt].ddc_pos;
    
    for ( ibpm = 1; ibpm < num_bpms; ibpm++ ) {

      if ( ana_cutfn( &(proc[ibpm][evt]) ) == BPM_GOOD_EVENT ) {

	residual -= proc[ibpm][evt].ddc_pos * coeffs[idx++];

	if (ANA_SVD_TILT & mode) 
	  residual -= proc[ibpm][evt].ddc_slope * coeffs[idx++];
      }
      else good_evt = 0;
    }

    if (!good_evt) continue;
    
    // Good event
    good_res[count] = residual - coeffs[idx];
    count++;
  }
  
  // Calculate the mean
  *mean = 0;
  for (evt = 0; evt < count; evt++) {
    *mean += good_res[evt];
  }

  *mean /= (double) count;

  // Calculate the RMS
  *rms = 0;
  for (evt = 0; evt < count; evt++) {

    *rms += (good_res[evt] - *mean) * (good_res[evt] - *mean);
  }

  *rms = sqrt( *rms / count );

  free(good_res);
  
  return BPM_SUCCESS;
}
