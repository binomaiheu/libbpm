/**
   @file
   @ingroup wave
*/

#include <bpm/bpm_wf.h>
#include <bpm/bpm_dsp.h>

complexwf_t* complexwf( int ns, double fs ) {
  
  complexwf_t *w;

  if ( fs <= 0 ) {
    bpm_error( "Cannot have sampling frequency <= 0. in complexwf()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  if ( ns > MAX_ALLOWED_NS ) {
    bpm_error( "Maximum allowed number of samples exceeded, failed to allocate.",
	       __FILE__, __LINE__ );
    return NULL;
  }

  if ( ns > 1 ) {
    w = (complexwf_t*) calloc( 1, sizeof(complexwf_t) );
    if ( ! w ) {
      bpm_error( "Cannot allocate memory for waveform structure in complexwf()",
		 __FILE__, __LINE__ );
      return NULL;
    }
    w->ns = ns;
    w->fs = fs;
    w->wf = (complex_t*) calloc( w->ns, sizeof( complex_t ) );
    if ( ! w->wf ) {
      bpm_error( "Cannot allocate memory for waveform data in complexwf()",
		 __FILE__, __LINE__ );
      free( w );
      return NULL;
    }
  } else {
    bpm_error( "Invalid number of samples in complexwf()", __FILE__, __LINE__ );
    return NULL;
  }

  return w;
}

// ----------------------------------------------------------------------------

complexwf_t* complexwf_copy_new( complexwf_t *w ) {
  
  int i = 0;
  complexwf_t *s = complexwf( w->ns, w->fs );
  
  if ( ! s ) {
    bpm_error( "Cannot allocate memory in complexwf_copy_new()", __FILE__, __LINE__ );
    return NULL;
  };

  for ( i=0; i<w->ns; i++ ) s->wf[i] = w->wf[i];

  return s;
}

// ----------------------------------------------------------------------------

int complexwf_copy( complexwf_t *copy, complexwf_t *src ) {
  
  int i = 0;

  if ( ! copy || ! src ) {
    bpm_error( "Invalid pointer arguments in complexwf_copy()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  };

  if ( complexwf_compat( copy, src ) ) {
    for ( i=0; i<copy->ns; i++ ) copy->wf[i] = src->wf[i];
  } else {
    bpm_error( "Incompatible waveforms for in complexwf_copy()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_subset( complexwf_t *sub, complexwf_t *w, int i1, int i2 ) {
  
  int i = 0;

  if ( ! sub || ! w ) {
    bpm_error( "Invalid pointer arguments in complexwf_subset()", __FILE__, __LINE__ );
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

int complexwf_setvalues( complexwf_t *w, complex_t *x ) {

  int i = 0;
  if ( ! w || ! x ) {
    bpm_error( "Invalid pointer arguments in complexwf_setvalues()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = x[i];

  return BPM_SUCCESS;
}
// ----------------------------------------------------------------------------

int complexwf_setfunction( complexwf_t *w, 
			   complex_t (*wffun)( double, int, double* ),
			   int npars, double *par ) {
  
  int i = 0;
  if ( ! w || ! wffun ) {
    bpm_error( "Invalid pointer arguments in complexwf_setfunction()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = (*wffun)( (double) i / w->fs, npars, par );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_reset( complexwf_t *w ) {

  int i = 0;

  if ( ! w  ) {
    bpm_error( "Invalid pointer argument in complexwf_reset()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = complex( 0., 0. );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

void complexwf_delete( complexwf_t *w ) {

  if ( w ) {
    if ( w->wf ) free( w->wf ); else
      bpm_warning( "Cannot free complexwf_t::wf pointer in complexwf()_delete, already NULL !",
		   __FILE__, __LINE__ );
    free( w );
  } else {
    bpm_warning( "Cannot free complexwf_t pointer in complexwf()_delete, already NULL !", 
		 __FILE__, __LINE__ );
  }
  
  return;
}

// ----------------------------------------------------------------------------

int complexwf_compat( complexwf_t *w1, complexwf_t *w2 ) {

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_compat()", 
	       __FILE__, __LINE__ );
    return 0;
  }
  
  return ((w1->ns==w2->ns)&&(fabs(w1->fs-w2->fs)<WF_EPS)?1:0);
}

// ----------------------------------------------------------------------------

int complexwf_add( complexwf_t *w1, complexwf_t *w2 ) {

  int i = 0;

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in complexwf_add()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! complexwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in complexwf_add()", __FILE__, __LINE__ );
  }
  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] = c_sum( w1->wf[i], w2->wf[i] );
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_subtract( complexwf_t *w1, complexwf_t *w2 ) {
  
  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in complexwf_subtract()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! complexwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in complexwf_subtract()", __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] = c_diff( w1->wf[i], w2->wf[i] );
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_multiply( complexwf_t *w1, complexwf_t *w2 ) {

  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in complexwf_multiply()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! complexwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in complexwf_multiply()", 
		 __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] = c_mult( w1->wf[i], w2->wf[i] );
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_divide( complexwf_t *w1, complexwf_t *w2 ) {
  
  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in complexwf_divide()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! complexwf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in complexwf_divide()", 
		 __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) {

    if ( c_isequal( w2->wf[i], complex( 0., 0. ) ) ) {
      bpm_warning( "Trapped division by 0+0i in complexwf_divide()", 
		   __FILE__, __LINE__ );
      w1->wf[i] = complex(0.,0.);
    } else {
      w1->wf[i] = c_div( w1->wf[i], w2->wf[i] );
    } 
    
  } 

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_scale( complex_t f, complexwf_t *w ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_scale()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = c_mult( f, w->wf[i] );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_bias( complex_t c, complexwf_t *w ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_bias()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = c_sum( c, w->wf[i] );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_add_cwtone( complexwf_t *w, double amp, double phase, double freq,
			 double phasenoise ) {
  
  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_add_cwtone()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    w->wf[i].re += amp * cos( 2. * PI * freq * (double) i /  w->fs + 
			      nr_rangauss( phase, phasenoise ) );
    
    w->wf[i].im += amp * sin( 2. * PI * freq * (double) i /  w->fs + 
			      nr_rangauss( phase, phasenoise ) );
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_add_dcywave( complexwf_t *w, double amp, double phase, double freq, 
			   double ttrig, double tdcy, double phasenoise ) {

  int    i = 0;
  double t = 0.;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_add_dcywave()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  for (i=0; i<w->ns; i++ ) {
    t = (double) i / w->fs;
    if ( t >= ttrig ) {
      w->wf[i].re += amp * exp( -( t - ttrig ) / tdcy ) * 
	cos( 2. * PI * freq * ( t - ttrig ) + nr_rangauss( phase, phasenoise ) );
      w->wf[i].im += amp * exp( -( t - ttrig ) / tdcy ) * 
	sin( 2. * PI * freq * ( t - ttrig ) + nr_rangauss( phase, phasenoise ) );
    }
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_add_noise( complexwf_t *w, double sigma ) {
  
  int i = 0;
  double amp, phase;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_add_noise()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    amp   = nr_rangauss( 0., sigma );
    phase = nr_ranuniform( 0., 2*PI );

    w->wf[i] = c_sum( w->wf[i], complex( amp*cos(phase), amp*sin(phase) ) );
  }
  
  return BPM_SUCCESS;

}

// ----------------------------------------------------------------------------

int complexwf_add_ampnoise( complexwf_t *w, double sigma ) {
  
  int i = 0;
  double amp, phase;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_add_ampnoise()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    amp   = nr_rangauss( c_abs( w->wf[i]), sigma );
    phase = c_arg( w->wf[i] );

    w->wf[i] = complex( amp * cos( phase ), amp * sin( phase ) );
  }

  return BPM_SUCCESS;

}

// ----------------------------------------------------------------------------

int complexwf_add_phasenoise( complexwf_t *w, double sigma ) {
  
  int i = 0;
  double amp, phase;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in complexwf_add_phasenoise()", 
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    amp      = c_abs( w->wf[i] ); 
    phase    = nr_rangauss( c_arg( w->wf[i]), sigma );

    w->wf[i] = complex( amp * cos( phase ), amp * sin( phase ) );
  }

  return BPM_SUCCESS;

}


// ----------------------------------------------------------------------------

void complexwf_print( FILE *of, complexwf_t *w ) {
  
  int i = 0;
  if ( !of || !w ) {
    bpm_error( "Invalid pointers in comlexwf_print()", __FILE__, __LINE__ );
    return;
  }
  
  fprintf( of, "Waveform:\n" );
  fprintf( of, "Number of samples  : %d\n",     w->ns );
  fprintf( of, "Sampling frequency : %f MHz\n", w->fs / _MHz__ );
  for( i = 0; i < w->ns; i++ )  fprintf( of, "  wf[%5d] = %.14e + i %.14e \n", 
					 i, w->wf[i].re, w->wf[i].im );
  fflush( of );

  return;
}

// ----------------------------------------------------------------------------

int complexwf_getreal( doublewf_t *re, complexwf_t *z ) {
  
  int i = 0;

  if ( ! re || ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getreal()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( re->ns != z->ns ) {
    bpm_warning( "Different number of samples in complex_getreal()",
		 __FILE__, __LINE__ );
  }
  
  for ( i=0; i < MIN(re->ns,z->ns); i++ ) re->wf[i] = z->wf[i].re;  

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_getimag( doublewf_t *im, complexwf_t *z ) {
  
  int i = 0;

  if ( ! im || ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getimag()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( im->ns != z->ns ) {
    bpm_warning( "Different number of samples in complex_getimag()",
		 __FILE__, __LINE__ );
  }
  
  for ( i=0; i < MIN(im->ns,z->ns); i++ ) im->wf[i] = z->wf[i].im;  

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_getamp( doublewf_t *r, complexwf_t *z ) {
  
  int i = 0;

  if ( ! r || ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getamp()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( r->ns != z->ns ) {
    bpm_warning( "Different number of samples in complex_getamp()",
		 __FILE__, __LINE__ );
  }
  
  for ( i=0; i < MIN(r->ns,z->ns); i++ ) r->wf[i] = c_abs( z->wf[i] );  

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_getphase( doublewf_t *theta, complexwf_t *z ) {
  
  int i = 0;

  if ( ! theta || ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getphase()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( theta->ns != z->ns ) {
    bpm_warning( "Different number of samples in complexwf_getphase()",
		 __FILE__, __LINE__ );
  }
  
  for ( i=0; i < MIN(theta->ns,z->ns); i++ ) {
    theta->wf[i] = c_arg( z->wf[i] );  
    norm_phase( &( theta->wf[i] ) );
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_setreal( complexwf_t *z, doublewf_t *re ){

  int i = 0;

  if ( ! re || ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_setreal()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( re->ns != z->ns ) {
    bpm_warning( "Different number of samples in complexwf_setreal()",
		 __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(re->ns,z->ns); i++ ) z->wf[i].re = re->wf[i];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int complexwf_setimag( complexwf_t *z, doublewf_t *im ){

  int i = 0;

  if ( ! im || ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_setreal()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  if ( im->ns != z->ns ) {
    bpm_warning( "Different number of samples in complexwf_setreal()",
		 __FILE__, __LINE__ );
  }
  
  for ( i=0; i < MIN(im->ns,z->ns); i++ ) z->wf[i].im = im->wf[i];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

doublewf_t* complexwf_getreal_new( complexwf_t *z ) {

  int i = 0;
  doublewf_t *w = NULL;

  if ( ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getreal_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  w = doublewf( z->ns, z->fs );
  if ( ! w ) {
    bpm_error( "Unable to allocate memory for waveform in complex_getreal_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  for (i=0;i<z->ns;i++) w->wf[i] = z->wf[i].re;

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* complexwf_getimag_new( complexwf_t *z ) {

  int i = 0;
  doublewf_t *w = NULL;

  if ( ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getimag_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  w = doublewf( z->ns, z->fs );
  if ( ! w ) {
    bpm_error( "Unable to allocate memory for waveform in complex_getimag_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  for (i=0;i<z->ns;i++) w->wf[i] = z->wf[i].im;

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* complexwf_getamp_new( complexwf_t *z ) {

  int i = 0;
  doublewf_t *w = NULL;

  if ( ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getamp_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  w = doublewf( z->ns, z->fs );
  if ( ! w ) {
    bpm_error( "Unable to allocate memory for waveform in complex_getamp_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  for (i=0;i<z->ns;i++) w->wf[i] = c_abs( z->wf[i] );

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* complexwf_getphase_new( complexwf_t *z ) {

  int i = 0;
  doublewf_t *w = NULL;

  if ( ! z ) {
    bpm_error( "Invalid pointer argument in complexwf_getphase_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  w = doublewf( z->ns, z->fs );
  if ( ! w ) {
    bpm_error( "Unable to allocate memory for waveform in complex_getphase_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  for (i=0;i<z->ns;i++) {
    w->wf[i] = c_arg( z->wf[i] );
    norm_phase( &(w->wf[i]) );
  }

  return w;
}


