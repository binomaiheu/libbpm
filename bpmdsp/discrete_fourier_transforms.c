/**
   @file
   @ingroup dsp
*/

#include "bpm/bpm_wf.h"
#include "bpm/bpm_dsp.h"

// prototype the two functions from fftsg.c that we're going to use 
// in this code.
void cdft(int, int, double *, int *, double *);
void rdft(int, int, double *, int *, double *);


// define the fft buffers
static int    *_fft_work_area        = NULL;  // bitreversal work area
static int     _fft_work_area_length = 0;     // ... buffer length 
static double *_fft_sc_table         = NULL;  // sin/cos value table
static int     _fft_sc_table_length  = 0;     // ... buffer length
static double *_fft_data             = NULL;  // fft input data array
static int     _fft_data_length      = 0;     // ... buffer length

// ----------------------------------------------------------------------------

int _is_pow2( int n ) {
  int p = 0, r = 0;
  do {
    r = n % 2;  n /= 2;  p++;
  } while ( r == 0 && n > 1 );
  if( r != 0 ) return FALSE;
  return p;
}

// ----------------------------------------------------------------------------

int _check_fft_buffers( int ns ) {

  // recommended work area buffer length
  int nw = 2+(1<<(int)(log(ns/2+0.5)/log(2))/2);
  int nt = ns/2;

  // check the work area buffer
  if ( _fft_work_area ) {
    if ( _fft_work_area_length < nw ) {
      bpm_warning( "FFT work buffer to small, increasing size...", __FILE__, __LINE__ );
      free( _fft_work_area );
      _fft_work_area = (int*) calloc( nw, sizeof(int) );
      if ( ! _fft_work_area ) {
	bpm_error( "Cannot allocate memory for FFT work buffer", __FILE__, __LINE__ );
	return BPM_FAILURE;
      }
      _fft_work_area_length = nw;
    }
  } else {
    bpm_warning( "Allocating FFT work buffer, no fft_initialise() found", __FILE__, __LINE__ );
    _fft_work_area = (int*) calloc( nw, sizeof(int) );
    if ( ! _fft_work_area ) {
      bpm_error( "Cannot allocate memory for FFT work buffer", __FILE__, __LINE__ );
      return BPM_FAILURE;
    }
    _fft_work_area_length = nw;
  }
  
  // check the sin/cos table buffer length
  if ( _fft_sc_table ) {
    if ( _fft_sc_table_length < nt ) {
      bpm_warning( "FFT sin/cos table too small, increasing size...", __FILE__, __LINE__ );
      free( _fft_sc_table );
      _fft_sc_table = (double*) calloc( nt, sizeof(double) );
      if ( ! _fft_sc_table ) {
	bpm_error( "Cannot allocate memory for FFT sin/cos table", __FILE__, __LINE__ );
	return BPM_FAILURE;
      }
      _fft_sc_table_length = nt;
    }
  } else {
    bpm_warning( "Allocating FFT sin/cos table buffer, no fft_initialise() found",
		 __FILE__, __LINE__ );
    _fft_sc_table = (double*) calloc( nt, sizeof(double) );
    if ( ! _fft_sc_table ) {
      bpm_error( "Cannot allocate memory for FFT sin/cos table", __FILE__, __LINE__ );
      return BPM_FAILURE;
    }
    _fft_sc_table_length = nt;
  }

  // check the fft data array buffer length
  if ( _fft_data ) {
    if ( _fft_data_length < 2*ns ) {
      bpm_warning( "FFT data buffer length too small, increasing size...", __FILE__, __LINE__ );
      free( _fft_data );
      _fft_data = (double*) calloc( 2*ns, sizeof(double) );
      if ( ! _fft_data ) {
	bpm_error( "Cannot allocate memory for FFT data buffer", __FILE__, __LINE__ );
	return BPM_FAILURE;
      }
      _fft_data_length = 2*ns;
    }
  } else {
    bpm_warning( "Allocating FFT data buffer, no fft_initialise() found",
		 __FILE__, __LINE__ );
    _fft_data = (double*) calloc( 2*ns, sizeof(double) );
    if ( ! _fft_data ) {
      bpm_error( "Cannot allocate memory for FFT data buffer", __FILE__, __LINE__ );
      return BPM_FAILURE;
    }
    _fft_data_length = 2*ns;
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------


int  fft_gen_tables( void ) {

  if ( _fft_work_area ) {
    _fft_work_area[0] = 0;
  } else {
    bpm_error( "FFT work buffer not allocated, cannot regenerate tables", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int fft_initialise( int ns ) {
  
  // recommended work area buffer length
  int nw = 2+(1<<(int)(log(ns/2+0.5)/log(2))/2);
  int nt = ns/2;

  if ( _fft_work_area || _fft_sc_table || _fft_data ) {
    bpm_error( "FFT buffers alread initialised, please cleanup first with fft_cleanup()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  _fft_work_area = (int*) calloc( nw, sizeof(int) );
  _fft_sc_table  = (double*) calloc( nt, sizeof(double) );
  _fft_data      = (double*) calloc( 2*ns, sizeof(double) );

  if ( ! _fft_work_area || ! _fft_sc_table || ! _fft_data ) {
    bpm_error( "Failed to allocate memory for the FFT buffers in fft_initialise",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  _fft_work_area_length = nw;
  _fft_sc_table_length  = nt;
  _fft_data_length      = 2*ns;
  // setting the first bit to 0 to tell the fft to 
  // generate the sin/cos table

  return fft_gen_tables();
}

// ----------------------------------------------------------------------------

void fft_cleanup( void ) {

  if( _fft_work_area ) free( _fft_work_area );
  if( _fft_sc_table )  free( _fft_sc_table );
  if( _fft_data )      free( _fft_data );

  _fft_work_area_length = 0;
  _fft_sc_table_length  = 0;
  _fft_data_length      = 0;

  return;
}

// ----------------------------------------------------------------------------

int complexfft( complexwf_t *z, int fft_mode ) {

  // in place fft on complex wf z 

  int i = 0;

  if (  ! z ) {
    bpm_error( "Invalid pointers in complexfft(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! _is_pow2(z->ns) ) {
    bpm_warning( "Number of samples is not of the form 2^n, may run into trouble with FFT !",
		 __FILE__, __LINE__ );
  }

  if ( _check_fft_buffers( z->ns ) == BPM_FAILURE ) {
    bpm_error( "Error checking FFT buffers in complexfft()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // fill FFT data buffer
  for ( i=0; i<z->ns; i++ ) {
    _fft_data[2*i]   = z->wf[i].re;
    _fft_data[2*i+1] = z->wf[i].im; 
  }

  // execute the cdft
  switch ( fft_mode ) {
    case FFT_FORWARD:
      cdft( 2*z->ns, 1, _fft_data, _fft_work_area, _fft_sc_table );
      break;
    case FFT_BACKWARD:
      cdft( 2*z->ns, -1, _fft_data, _fft_work_area, _fft_sc_table );
      break;
    default:
      bpm_error( "Unknown FFT mode in complexfft()", __FILE__, __LINE__ );
      return BPM_FAILURE;
      break;
  }
  
  // reset the data
  for ( i=0; i<z->ns; i++ ) {
    z->wf[i].re = _fft_data[2*i];
    z->wf[i].im = _fft_data[2*i+1];
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int realfft( doublewf_t *y, int fft_mode, complexwf_t *z ) {

  // FFT_FORWARD : real wf of ns samples y -> complex wf of ns/2 samples z
  // FFT_BACKWARD: complex wf of ns/2 samples in z -> real wf of ns samples in y

  int i = 0;

  if ( ! y || ! z ) {
    bpm_error( "Invalid pointers in realfft(...)", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! _is_pow2(z->ns) ) {
    bpm_warning( "Number of samples is not of the form 2^n, may run into trouble with FFT !",
		 __FILE__, __LINE__ );
  }

  if ( _check_fft_buffers( z->ns ) == BPM_FAILURE ) {
    bpm_error( "Error checking FFT buffers in complexfft()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }


  // execute the cdft
  switch ( fft_mode ) {
    case FFT_FORWARD:
      // fill FFT data buffer
      for ( i=0; i<z->ns; i++ ) _fft_data[i] = y->wf[i];

      rdft( z->ns, 1, _fft_data, _fft_work_area, _fft_sc_table );

      // refill the data
      for ( i=0; i<z->ns/2; i++ ) {
	z->wf[i].re = z->wf[z->ns-1-i].re = _fft_data[2*i];
	z->wf[i].im = z->wf[z->ns-1-i].im = _fft_data[2*i+1];
      }

      break;
    case FFT_BACKWARD:
      // fill FFT data buffer
      for ( i=0; i<z->ns/2; i++ ) {
	_fft_data[2*i]   = z->wf[i].re;
	_fft_data[2*i+1] = z->wf[i].im;
      }

      rdft( z->ns, -1, _fft_data, _fft_work_area, _fft_sc_table );
      
      // refill the data
      for ( i=0; i<z->ns; i++ ) y->wf[i] = _fft_data[i];
      break;

    default:
      bpm_error( "Unknown FFT mode in complexfft()", __FILE__, __LINE__ );
      return BPM_FAILURE;
      break;
  }
  

  return BPM_SUCCESS;
}
