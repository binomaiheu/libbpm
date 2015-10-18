/**
   @file 
   @ingroup dsp
*/


#include "bpm/bpm_dsp.h"

// the ddc buffers
static doublewf_t *_ddc_buffer_real = NULL;
static doublewf_t *_ddc_buffer_imag = NULL;

int _check_ddc_buffers( int ns, double fs ) {

  if ( _ddc_buffer_real ) {
    if ( ( _ddc_buffer_real->ns != ns ) || ( ( _ddc_buffer_real->fs - fs ) > 1.0e-10 ) ) {
      bpm_warning( "Reallocating _ddc_buffer_real with different number of samples & fs!",
		   __FILE__, __LINE__ );
      doublewf_delete ( _ddc_buffer_real );
      _ddc_buffer_real = doublewf( ns, fs );
    }
  } else {
    bpm_warning( "Allocating DDC-Re buffer, no ddc_initialise() found", __FILE__, __LINE__ );
    _ddc_buffer_real = doublewf( ns, fs );
  }

  if ( _ddc_buffer_imag ) {
    if ( ( _ddc_buffer_imag->ns != ns ) || ( ( _ddc_buffer_imag->fs - fs ) > 1.0e-10 ) ) {
      bpm_warning( "Reallocating _ddc_buffer_imag with different number of samples & fs!",
		   __FILE__, __LINE__ );
      doublewf_delete ( _ddc_buffer_imag );
      _ddc_buffer_imag = doublewf( ns, fs );
    }
  } else {
    bpm_warning( "Allocating DDC-Im  buffer, no ddc_initialise() found", __FILE__, __LINE__ );
    _ddc_buffer_imag = doublewf( ns, fs );
  }

  if ( ! _ddc_buffer_real || ! _ddc_buffer_imag ) {
    bpm_error( "Cannot (re-)allocate memory for DDC buffers :(!", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int ddc_initialise( int ns, double fs ) {

  if ( _ddc_buffer_real || _ddc_buffer_imag ) {
    bpm_error( "DDC buffers already existing, please cleanup first with ddc_cleanup() !", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  _ddc_buffer_real = doublewf( ns, fs );
  _ddc_buffer_imag = doublewf( ns, fs );
  if( ! _ddc_buffer_real || ! _ddc_buffer_imag ) {
    bpm_error( "Failed to allocate memory for DDC buffers", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

void ddc_cleanup( void ) {
  if ( _ddc_buffer_real ) doublewf_delete( _ddc_buffer_real );
  if ( _ddc_buffer_imag ) doublewf_delete( _ddc_buffer_imag );
  return;
}

// ----------------------------------------------------------------------------

int ddc( doublewf_t *w, double f, filter_t *filter, complexwf_t *dcw,
	 doublewf_t *bufre,      // Optional buffer for real part
	 doublewf_t *bufim ) {   // Optional buffer for imag part
  
  int i = 0;
  doublewf_t *re, *im; // these are the buffers we'll use, either user or internal !

  if ( ! ( bufre && bufim ) ) {
    if ( _check_ddc_buffers( dcw->ns, dcw->fs ) ) return BPM_FAILURE;
    re = _ddc_buffer_real;
    im = _ddc_buffer_imag;
  } else {
    re = bufre;
    im = bufim;
  }
  
  // downmix entire waveform...
  for (i=0;i<w->ns;i++){
    re->wf[i] = w->wf[i] * cos( 2.*PI*f * (double) i / w->fs );
    im->wf[i] = w->wf[i] * sin( 2.*PI*f * (double) i / w->fs );
  }
  
  // filter 2 omega component
  if( apply_filter( filter, re ) == BPM_FAILURE ) return BPM_FAILURE;
  if( apply_filter( filter, im ) == BPM_FAILURE ) return BPM_FAILURE;
  
  // return complex waveform
  complexwf_setreal( dcw, re );
  complexwf_setimag( dcw, im );
  
  return BPM_SUCCESS;
}
