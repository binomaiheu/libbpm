/**
   @file
   @ingroup wave
*/

#include <bpm/bpm_wf.h>

int wfstat_reset( wfstat_t *s ) {

  if ( ! s ) {
    bpm_error( "Invalid pointer argument in reset_wfstats()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // re-initialise the structure
  s->imax = 0;
  s->imin = 0;
  s->max  = -DBL_MAX;
  s->min  =  DBL_MAX;
  s->mean = 0.;
  s->rms  = 0.;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

void wfstat_print( FILE *of, wfstat_t *s ) {

  if ( ! of || ! s ) {
    bpm_error( "Invalid pointer arguments in wfstat_print()", __FILE__, __LINE__ );
    return;
  }

  fprintf( of, "Basic waveform statistics:\n" );
  fprintf( of, " - maximum ... : wf[%d] = %.14e\n", s->imax, s->max );
  fprintf( of, " - minimum ... : wf[%d] = %.14e\n", s->imin, s->min );
  fprintf( of, " - mean ...... : %.14e\n", s->mean );
  fprintf( of, " - rms ....... : %.14e\n", s->rms );
  fflush( of );

  return;
}
