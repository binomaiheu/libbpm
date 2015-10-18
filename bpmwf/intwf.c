/**
   @file
   @ingroup wave
*/

#include <bpm/bpm_wf.h>

intwf_t* intwf( int ns, double fs ) {
  
  intwf_t *w;

  if ( fs <= 0 ) {
    bpm_error( "Cannot have sampling frequency <= 0. in intwf()",
	       __FILE__, __LINE__ );
    return NULL;
  }
  
  if ( ns > MAX_ALLOWED_NS ) {
    bpm_error( "Maximum allowed number of samples exceeded, failed to allocate.",
               __FILE__, __LINE__ );
    return NULL;
  }

 if ( ns > 1 ) {
    w = (intwf_t*) calloc( 1, sizeof(intwf_t) );
    if ( ! w ) {
      bpm_error( "Cannot allocate memory for waveform structure in intwf()",
		 __FILE__, __LINE__ );
      return NULL;
    }
    w->ns = ns;
    w->fs = fs;
    w->wf = (int*) calloc( w->ns, sizeof( int ) );
    if ( ! w->wf ) {
      bpm_error( "Cannot allocate memory for waveform data in intwf()",
		 __FILE__, __LINE__ );
      free( w );
      return NULL;
    }
  } else {
    bpm_error( "Invalid number of samples in intwf()", __FILE__, __LINE__ );
    return NULL;
  }

  return w;
}

// ----------------------------------------------------------------------------

intwf_t* intwf_sample_series( int ns, double fs ) {

  int i = 0;
  intwf_t *w = intwf( ns, fs );

  if ( ! w ) return w;
  for ( i=0; i<w->ns; i++ ) w->wf[i] = i;

  return w;
}

// ----------------------------------------------------------------------------

intwf_t* intwf_copy_new( intwf_t *w ) {
  
  int i = 0;
  intwf_t *s = intwf( w->ns, w->fs );
  
  if ( ! s ) {
    bpm_error( "Cannot allocate memory in intwf_copy_new()", 
	       __FILE__, __LINE__ );
    return NULL;
  };

  for ( i=0; i<w->ns; i++ ) s->wf[i] = w->wf[i];

  return s;
}

// ----------------------------------------------------------------------------

int intwf_copy( intwf_t *copy, intwf_t *src ) {
  
  int i = 0;

  if ( ! copy || ! src ) {
    bpm_error( "Invalid pointer arguments in intwf_copy()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  };

  if ( intwf_compat( copy, src ) ) {
    for ( i=0; i<copy->ns; i++ ) copy->wf[i] = src->wf[i];
  } else {
    bpm_error( "Incompatible waveforms for in intwf_copy()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_subset( intwf_t *sub, intwf_t *w, int i1, int i2 ) {
  
  int i = 0;

  if ( ! sub || ! w ) {
    bpm_error( "Invalid pointer arguments in intwf_subset()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  };

  // reset the number of samples in the subset
  sub->ns = 0;
  sub->fs = w->fs;

  for ( i=MAX(0,i1); i<=MIN(w->ns-1,i2); i++ ) {
    sub->wf[i] = w->wf[i-i1];
    sub->ns++;
  }


  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_setvalues( intwf_t *w, int *x ) {

  int i = 0;
  if ( ! w || ! x ) {
    bpm_error( "Invalid pointer arguments in intwf_setvalues()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = x[i];

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_setfunction( intwf_t *w, 
		       int (*wffun)( double, int, double* ),
		       int npars, double *par ) {
  
  int i = 0;
  if ( ! w || ! wffun ) {
    bpm_error( "Invalid pointer arguments in intwf_setfunction()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  for ( i=0; i<w->ns; i++ ) w->wf[i] = (*wffun)( (double) i / w->fs, npars, par );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_reset( intwf_t *w ) {

  int i = 0;

  if ( ! w  ) {
    bpm_error( "Invalid pointer argument in intwf_reset()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = 0;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

void intwf_delete( intwf_t *w ) {

  if ( w ) {
    if ( w->wf ) free( w->wf ); else
      bpm_warning( "Cannot free intwf_t::wf pointer in intwf()_delete, already NULL !",
                   __FILE__, __LINE__ );
    free( w );
  } else {
    bpm_warning( "Cannot free intwf_t pointer in intwf()_delete, already NULL !", 
                 __FILE__, __LINE__ );
  }
  
  return;
}

// ----------------------------------------------------------------------------

doublewf_t* doublewf_cast_new( intwf_t *iw ) {

  int i = 0;
  doublewf_t *w;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_cast_new()", 
               __FILE__, __LINE__ );
    return NULL;
  }

  w = doublewf( iw->ns, iw->fs );
  if ( ! w ) {
    bpm_error( "Cannot allocate memory for doublewf_t in doublewf_cast_new()",
               __FILE__, __LINE__ );
    return NULL;
  }
  
  // cast to the double waveform
  for (i=0;i<iw->ns;i++) w->wf[i] = (double) iw->wf[i];

  return w;
}

// ----------------------------------------------------------------------------

int doublewf_cast( doublewf_t *w, intwf_t *iw ) {

  int i = 0;

  if ( ! w || ! iw ) {
    bpm_error( "Invalid pointer argument in doublewf_cast()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // cast to the double waveform
  for (i=0;i<iw->ns;i++) w->wf[i] = (double) iw->wf[i];

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_compat( intwf_t *w1, intwf_t *w2 ) {

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in intwf_compat()", 
               __FILE__, __LINE__ );
    return 0;
  }
  
  return ((w1->ns==w2->ns)&&(fabs(w1->fs-w2->fs)<WF_EPS)?1:0);
}

// ----------------------------------------------------------------------------

int intwf_add( intwf_t *w1, intwf_t *w2 ) {

  int i = 0;

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in intwf_add()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! intwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in intwf_add()", __FILE__, __LINE__ );
  }
  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] += w2->wf[i];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_subtract( intwf_t *w1, intwf_t *w2 ) {
  
  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in intwf_subtract()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! intwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in intwf_subtract()", __FILE__, __LINE__ );
  }
  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] -= w2->wf[i];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_multiply( intwf_t *w1, intwf_t *w2 ) {

  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in intwf_multiply()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! intwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in intwf_multiply()", 
                 __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] *= w2->wf[i];

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_divide( intwf_t *w1, intwf_t *w2 ) {
  
  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in intwf_divide()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! intwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in intwf_divide()", 
                 __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) {
    if ( w2->wf[i] != 0. ) {
      w1->wf[i] /= w2->wf[i] ;
    } else {
      bpm_warning( "Trapped division by 0. in intwf_divide()", 
                   __FILE__, __LINE__ );
      w1->wf[i] = 0.;
    } 
  } 

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_scale( int f, intwf_t *w ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_scale()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] *= f;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_bias( int c, intwf_t *w ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_bias()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] += c;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_add_cwtone( intwf_t *w, double amp, double phase, double freq,
		      double phasenoise ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_add_cwtone()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  for (i=0; i<w->ns; i++ ) 
    w->wf[i] += (int) dround ( amp * cos( 2. * PI * freq * (double) i /  w->fs + 
					  nr_rangauss( phase, phasenoise ) ) );
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_add_dcywave( intwf_t *w, double amp, double phase, double freq, 
		       double ttrig, double tdcy, double phasenoise ) {
  
  int    i = 0;
  double t = 0.;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_add_dcywave()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  for (i=0; i<w->ns; i++ ) {
    t = (double) i / w->fs;
    if ( t >= ttrig ) {
      w->wf[i] += (int) dround( amp * exp( -( t - ttrig ) / tdcy ) * 
				cos( 2. * PI * freq * ( t - ttrig ) + 
				     nr_rangauss( phase, phasenoise ) ) );
    }
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_add_ampnoise( intwf_t *w, double sigma ) {
  
  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_add_ampnoise()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    w->wf[i] += (int) dround( nr_rangauss( 0., sigma ) );
  }

  return BPM_SUCCESS;

}

// ----------------------------------------------------------------------------

int intwf_basic_stats( intwf_t *w, int s0, int s1, wfstat_t *stats ) {
  
  doublewf_t *dw = NULL;
  
  if ( ! w || ! stats ) {
    bpm_error( "Invalid pointer arguments in intwf_basic_stats()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  dw = doublewf_cast_new( w );
  if ( ! dw ) {
    bpm_error( "Cannot allocate memory for temporary doublewf in intwf_basic_stats",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( doublewf_basic_stats( dw, s0, s1, stats ) ) return BPM_FAILURE;
  
  doublewf_delete( dw );
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_derive( intwf_t *w ) {

  int     i = 0;
  double dt = 0;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_derive()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // the time step...
  dt = 1./w->fs;

  // derivative...
  for ( i=0; i<w->ns-1; i++ ) w->wf[i] = (int) dround( ( w->wf[i+1] - w->wf[i] ) / dt );
  
  // the last sample.. extrapolate linear, so same slope
  w->wf[w->ns-1] = w->wf[w->ns-2];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int intwf_integrate( intwf_t *w ) {

  double tmp = 0., tmp_prev = 0.;
  double dt = 0.;
  int i = 0;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_integrate()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // the time step...
  dt = 1./w->fs;

  // integral
  for ( i=0; i<w->ns; i++ ) {
    tmp  = (double) w->wf[i];
    tmp *= dt;
    
    if ( i > 0 ) tmp += tmp_prev;

    w->wf[i] = (int) dround( tmp );
    tmp_prev = tmp;
  }
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

 void intwf_print( FILE *of, intwf_t *w ) {

  int i = 0;
  if ( !of || !w ) {
    bpm_error( "Invalid pointers in intwf_print()", __FILE__, __LINE__ );
    return;
  }
  
  fprintf( of, "Waveform:\n" );
  fprintf( of, "Number of samples  : %d\n",     w->ns );
  fprintf( of, "Sampling frequency : %f MHz\n", w->fs / _MHz__ );
  for( i = 0; i < w->ns; i++ )  fprintf( of, "  wf[%5d] = %d \n", i, w->wf[i] );
  fflush( of );

  return;
}

// ----------------------------------------------------------------------------

int intwf_getvalue( intwf_t *w, double t, unsigned int mode ) {
  
  int val = 0;
  doublewf_t *dw = NULL;

  if ( ! w ) {
    bpm_error( "Invalid pointer arguments in intwf_getvalue()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  dw = doublewf_cast_new( w );
  if ( ! dw ) {
    bpm_error( "Cannot allocate memory for temporary doublewf in intwf_getvalue()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  val = (int) dround ( doublewf_getvalue( dw, t, mode ) );

  doublewf_delete( dw );
  
  return val;
}

// ----------------------------------------------------------------------------

int intwf_resample( intwf_t *w2, double fs, intwf_t *w1, unsigned int mode ) {

  int i = 0;
  doublewf_t *dw = NULL;

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in intwf_resample()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  dw = doublewf_cast_new( w1 );
  if ( ! dw ) {
    bpm_error( "Cannot allocate memory for temporary doublewf in intwf_resample()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // number of samples
  w2->ns = (int) (w1->ns * fs / w1->fs);
  w2->fs = fs;

  if ( w2->ns > MAX_ALLOWED_NS ) {
    bpm_error( "Maximum allowed number of samples exceeded in intwf_resample()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( w2->ns < 1 ) {
    bpm_error( "Number of new samples is zero in intwf_resample()",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // fill the new waveform
  for ( i=0; i<w2->ns; i++ ) {
    w2->wf[i] = (int) dround( doublewf_getvalue( dw, (double) i / w2->fs, mode ) );
  }

  doublewf_delete( dw );

  return BPM_SUCCESS;
}
