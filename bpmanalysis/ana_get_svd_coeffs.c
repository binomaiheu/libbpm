/**
   @file
   @ingroup analysis
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_analysis.h>
#include <bpm/bpm_nr.h>

int ana_get_svd_coeffs( bpmproc_t **proc, int num_bpms, int num_svd, 
			int total_num_evts, double *coeffs, int mode ) {
  

  //----------------------------------------------------------
  // Set up mode specfic info
  int num_pars = num_bpms;  // At least the constant term + num_bpms

  if (ANA_SVD_TILT & mode) num_pars += num_bpms - 1;

  //----------------------------------------------------------
  // Loop over the procs, apply the cuts and store the data in a gsl_matrix
  int evt, ibpm, idx, count = 0, good_evt;

  // Set up matrices
  gsl_matrix *A = gsl_matrix_calloc( num_svd, num_pars );
  gsl_matrix *V = gsl_matrix_calloc(num_pars, num_pars);
  gsl_vector *S = gsl_vector_calloc(num_pars);
  gsl_vector *work = gsl_vector_calloc(num_pars);
  gsl_vector *b = gsl_vector_calloc(num_svd);
  gsl_vector *x = gsl_vector_calloc(num_pars);

  for ( evt = 0; evt < total_num_evts; evt++ ) {

    idx = 0;
    good_evt = 1;

    for ( ibpm = 1; ibpm < num_bpms; ibpm++ ) {

      if ( ana_cutfn( &(proc[ibpm][evt]) ) == BPM_GOOD_EVENT ) {
      		
	gsl_matrix_set( A, count, idx++, proc[ibpm][evt].ddc_pos );

	if (ANA_SVD_TILT & mode) 
	  gsl_matrix_set( A, count, idx++, proc[ibpm][evt].ddc_slope );
      }
      else good_evt = 0;
    }

    if (!good_evt) continue;
    
    // Good event
    gsl_matrix_set( A, count, idx, 1 );
    gsl_vector_set( b, count, proc[0][evt].ddc_pos );
    count++;
  }

  // if count != evt

  // Perform the SVD
  gsl_linalg_SV_decomp ( A, V, S, work);
  gsl_linalg_SV_solve (A, V, S, b, x);

  // Output the coeffs
  for (idx = 0; idx < num_pars; idx++) {
    coeffs[idx] = gsl_vector_get(x, idx);
  }
  
  return BPM_SUCCESS;
}
