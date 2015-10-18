/**
   @file
   @ingroup wave
*/

#include  <bpm/bpm_wf.h>

doublewf_t* doublewf( int ns, double fs ) {
  
  doublewf_t *w;

  if ( fs <= 0 ) {
    bpm_error( "Cannot have sampling frequency <= 0. in doublewf()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  if ( ns > MAX_ALLOWED_NS ) {
    bpm_error( "Maximum allowed number of samples exceeded, failed to allocate.",
	       __FILE__, __LINE__ );
    return NULL;
  }

  if ( ns > 1 ) {
    w = (doublewf_t*) calloc( 1, sizeof(doublewf_t) );
    if ( ! w ) {
      bpm_error( "Cannot allocate memory for waveform structure in doublewf()",
		 __FILE__, __LINE__ );
      return NULL;
    }
    w->ns = ns;
    w->fs = fs;
    w->wf = (double*) calloc( w->ns, sizeof( double ) );
    if ( ! w->wf ) {
      bpm_error( "Cannot allocate memory for waveform data in doublewf()",
		 __FILE__, __LINE__ );
      free( w );
      return NULL;
    }
  } else {
    bpm_error( "Invalid number of samples in doublewf()", __FILE__, __LINE__ );
    return NULL;
  }

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* doublewf_sample_series( int ns, double fs ) {

  int i = 0;
  doublewf_t *w = doublewf( ns, fs );

  if ( ! w ) return w;
  for ( i=0; i<w->ns; i++ ) w->wf[i] = (double) i;

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* doublewf_time_series( int ns, double fs ) {

  int i = 0;
  doublewf_t *w = doublewf( ns, fs );

  if ( ! w ) return w;
  for ( i=0; i<w->ns; i++ ) w->wf[i] = (double) i / w->fs;

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* doublewf_frequency_series( int ns, double fs ) {

  int i = 0;
  doublewf_t *w = doublewf( ns, fs );

  if ( ! w ) return w;
  for ( i=0; i<w->ns; i++ ) w->wf[i] = (double) i * w->fs / (double) w->ns;

  return w;
}

// ----------------------------------------------------------------------------

doublewf_t* doublewf_copy_new( doublewf_t *w ) {
  
  int i = 0;
  doublewf_t *s = doublewf( w->ns, w->fs );
  
  if ( ! s ) {
    bpm_error( "Cannot allocate memory in doublewf_copy_new()", __FILE__, __LINE__ );
    return NULL;
  };

  for ( i=0; i<w->ns; i++ ) s->wf[i] = w->wf[i];

  return s;
}

// ----------------------------------------------------------------------------

int doublewf_copy( doublewf_t *copy, doublewf_t *src ) {
  
  int i = 0;

  if ( ! copy || ! src ) {
    bpm_error( "Invalid pointer arguments in doublewf_copy()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  };

  if ( doublewf_compat( copy, src ) ) {
    for ( i=0; i<copy->ns; i++ ) copy->wf[i] = src->wf[i];
  } else {
    bpm_error( "Incompatible waveforms for in doublewf_copy()", __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_subset( doublewf_t *sub, doublewf_t *w, int i1, int i2 ) {
  
  int i = 0;

  if ( ! sub || ! w ) {
    bpm_error( "Invalid pointer arguments in doublewf_subset()", __FILE__, __LINE__ );
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

int doublewf_setvalues( doublewf_t *w, double *x ) {

  int i = 0;
  if ( ! w || ! x ) {
    bpm_error( "Invalid pointer arguments in doublewf_setvalues()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = x[i];

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_setfunction( doublewf_t *w, 
			  double (*wffun)( double, int, double* ),
			  int npars, double *par ) {

  int i = 0;
  if ( ! w || ! wffun ) {
    bpm_error( "Invalid pointer arguments in doublewf_setfunction()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = (*wffun)( (double) i / w->fs, npars, par );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_reset( doublewf_t *w ) {

  int i = 0;

  if ( ! w  ) {
    bpm_error( "Invalid pointer argument in doublewf_reset()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] = 0.;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

void doublewf_delete( doublewf_t *w ) {

  if ( w ) {
    if ( w->wf ) free( w->wf ); else
      bpm_warning( "Cannot free doublewf_t::wf pointer in doublewf()_delete, already NULL !",
		   __FILE__, __LINE__ );
    free( w );
  } else {
    bpm_warning( "Cannot free doublewf_t pointer in doublewf()_delete, already NULL !", 
		 __FILE__, __LINE__ );
  }
  
  return;
}

// ----------------------------------------------------------------------------

intwf_t* intwf_cast_new( doublewf_t *w ) {

  int      i = 0;
  intwf_t *iw;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in intwf_cast_new()", 
	       __FILE__, __LINE__ );
    return NULL;
  }

  iw = intwf( w->ns, w->fs );
  if ( ! iw ) {
    bpm_error( "Cannot allocate memory for intwf_t in intwf_cast_new()",
	       __FILE__, __LINE__ );
    return NULL;
  }

  // cast to the integer waveform
  for (i=0;i<iw->ns;i++) iw->wf[i] = (int) dround( w->wf[i] );

  return iw;
}

// ----------------------------------------------------------------------------

int intwf_cast( intwf_t *iw, doublewf_t *w ) {

  int      i = 0;

  if ( ! w || ! iw ) {
    bpm_error( "Invalid pointer argument in intwf_cast()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // cast to the integer waveform
  for (i=0;i<iw->ns;i++) iw->wf[i] = (int) dround( w->wf[i] );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_compat( doublewf_t *w1, doublewf_t *w2 ) {

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_compat()", 
	       __FILE__, __LINE__ );
    return 0;
  }
  
  return ((w1->ns==w2->ns)&&(fabs(w1->fs-w2->fs)<WF_EPS)?1:0);
}

// ----------------------------------------------------------------------------

int doublewf_add( doublewf_t *w1, doublewf_t *w2 ) {

  int i = 0;

  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_add()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! doublewf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in doublewf_add()", __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] += w2->wf[i];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_subtract( doublewf_t *w1, doublewf_t *w2 ) {
  
  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_subtract()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! doublewf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in doublewf_subtract()", __FILE__, __LINE__ );
  }
  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] -= w2->wf[i];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_multiply( doublewf_t *w1, doublewf_t *w2 ) {

  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_multiply()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! doublewf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in doublewf_multiply()", 
		 __FILE__, __LINE__ );
  }
  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) w1->wf[i] *= w2->wf[i];

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_divide( doublewf_t *w1, doublewf_t *w2 ) {
  
  int i = 0;
  
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_divide()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( ! doublewf_compat( w1, w2 ) ) {
    bpm_warning( "Incompatible waveforms in doublewf_divide()", 
		 __FILE__, __LINE__ );
  }

  for ( i=0; i < MIN(w1->ns,w2->ns); i++ ) {
    if ( w2->wf[i] != 0. ) {
      w1->wf[i] /= w2->wf[i] ;
    } else {
      bpm_warning( "Trapped division by 0. in doublewf_divide()", 
		   __FILE__, __LINE__ );
      w1->wf[i] = 0.;
    } 
  } 

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_scale( double f, doublewf_t *w ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_scale()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] *= f;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_bias( double c, doublewf_t *w ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_bias()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for ( i=0; i<w->ns; i++ ) w->wf[i] += c;

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_add_cwtone( doublewf_t *w, double amp, double phase, double freq,
			 double phasenoise ) {

  int i = 0;
  
  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_add_cwtone()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) 
    w->wf[i] += amp * cos( 2. * PI * freq * (double) i /  w->fs + 
			   nr_rangauss( phase, phasenoise ) );
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_add_dcywave( doublewf_t *w, double amp, double phase, double freq, 
			  double ttrig, double tdcy, double phasenoise ) {

  int    i = 0;
  double t = 0.;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_add_dcywave()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    t = (double) i / w->fs;
    if ( t >= ttrig ) {
      w->wf[i] += amp * exp( -( t - ttrig ) / tdcy ) * 
	cos( 2. * PI * freq * ( t - ttrig ) + nr_rangauss( phase, phasenoise ) );
    }
  }

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_add_ampnoise( doublewf_t *w, double sigma ) {
  
  int i = 0;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_add_ampnoise()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  for (i=0; i<w->ns; i++ ) {
    w->wf[i] += nr_rangauss( 0., sigma );
  }

  return BPM_SUCCESS;

}

// ----------------------------------------------------------------------------

int doublewf_basic_stats( doublewf_t *w, int s0, int s1, wfstat_t *stats ) {
  
  int i = 0;

  if ( ! w || ! stats ) {
    bpm_error( "Invalid pointer arguments in doublewf_basic_stats()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // reset statistics
  wfstat_reset( stats );

  // check the limits
  if ( s0 > s1 ) {
    bpm_warning( "Swapping limits in doublewf_basic_stats()", __FILE__, __LINE__ );
    i=s1; s1=s0; s0=i;
  }

  // check boundaries
  if ( s0 <  0     ) s0 = 0;
  if ( s1 >= w->ns ) s1 = w->ns-1;

  for (i=s0; i<=s1; i++ ) {

    stats->mean += w->wf[i];
    stats->rms  += SQR( w->wf[i] );

    if ( w->wf[i] > stats->max ) { stats->max = w->wf[i]; stats->imax = i; }
    if ( w->wf[i] < stats->min ) { stats->min = w->wf[i]; stats->imin = i; }
  }

  stats->mean /= (double) ( s1 - s0 + 1.);
  stats->rms   = sqrt( stats->rms / (double) ( s1 - s0 + 1.) - SQR( stats->mean ) );

  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_derive( doublewf_t *w ) {

  int     i = 0;
  double dt = 0;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_derive()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }
  
  // the time step...
  dt = 1./w->fs;

  // derivative...
  for ( i=0; i<w->ns-1; i++ ) w->wf[i] = ( w->wf[i+1] - w->wf[i] ) / dt;
  
  // the last sample.. extrapolate linear, so same slope
  w->wf[w->ns-1] = w->wf[w->ns-2];
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------

int doublewf_integrate( doublewf_t *w ) {

  int     i = 0;
  double dt = 0;

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_integrate()", 
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // the time step...
  dt = 1./w->fs;

  // integral
  for ( i=0; i<w->ns; i++ ) {
    w->wf[i] *= dt;
    if ( i > 0 ) w->wf[i] += w->wf[i-1];
  }
  
  return BPM_SUCCESS;
}

// ----------------------------------------------------------------------------
void doublewf_print( FILE *of, doublewf_t *w ) {

  int i = 0;
  if ( !of || !w ) {
    bpm_error( "Invalid pointers in doublewf_print()", __FILE__, __LINE__ );
    return;
  }

  fprintf( of, "Waveform:\n" );
  fprintf( of, "Number of samples  : %d\n",     w->ns );
  fprintf( of, "Sampling frequency : %f MHz\n", w->fs / _MHz__ );
  for( i = 0; i < w->ns; i++ )  fprintf( of, "  wf[%5d] = %.14e \n", i, w->wf[i] );
  fflush( of );

  return;
}

// ----------------------------------------------------------------------------

double doublewf_getvalue( doublewf_t *w, double t, unsigned int mode ) {
  
  double val = 0.;
  int i = 0, i0 = 0, i1 = 0;
  
  int lanczos_a = 3; // order or something, kernel width...

  if ( ! w ) {
    bpm_error( "Invalid pointer argument in doublewf_sample()", 
	       __FILE__, __LINE__ );
    return -DBL_MAX;
  }
  
  // lanczos(x) kernel interpolate, best for waveform reconstrucution
  if ( mode & WF_LANCZOS ) {
    for ( i=0; i<w->ns; i++ ) {
      val += w->wf[i] * lanczos( ( t - (double) i / w->fs ) * w->fs, lanczos_a );
    }
    return val;

  } 

  // sinc(x) kernel interpolate, good for waveform reconstrucution
  if ( mode & WF_SINC ) {
    for ( i=0; i<w->ns; i++ ) {
      val += w->wf[i] * sinc( ( t - (double) i / w->fs ) * w->fs );
    }
    return val;

  } 

  // determine the two closest samples
  i0 = (int) ( t * w->fs );
  i1 = i0 + 1;
  
  if ( i0 < 0 ) i0 = 0;
  if ( i0 > w->ns-2 ) i0 = w->ns-2;

  if ( i1 < 1 ) i1 = 1;
  if ( i1 > w->ns-1 ) i1 = w->ns-1;
  

  // linear interpolation
  if ( mode & WF_LINEAR ) {
    val = w->wf[i0] + ( w->wf[i1] - w->wf[i0] ) * ( t * w->fs - (double) i0 );
    return val;
  }


  // parabolic interpolation
  if ( mode & WF_QUADRATIC ) {

    if ( ( t*w->fs - (double) i0 ) < 0.5 ) {
      if ( i0 > 0 ) {
	return nr_quadinterpol( t, 
				(double)(i0-1)/w->fs, (double)i0/w->fs, (double)(i1)/w->fs,
				w->wf[i0-1], w->wf[i0], w->wf[i1] );
      } else {
	return nr_quadinterpol( t, 
				(double)i0/w->fs, (double)i1/w->fs, (double)(i1+1)/w->fs,
				w->wf[i0], w->wf[i1], w->wf[i1+1] );
      }

    } else {

      if ( i1 < w->ns-1 ) {
	return nr_quadinterpol( t, 
				(double)i0/w->fs, (double)i1/w->fs, (double)(i1+1)/w->fs,
				w->wf[i0], w->wf[i1], w->wf[i1+1] );
      } else {
	return nr_quadinterpol( t, 
				(double)(i0-1)/w->fs, (double)i0/w->fs, (double)(i1)/w->fs,
				w->wf[i0-1], w->wf[i0], w->wf[i1] );
      }

    }

  }

  // no interpolation ( WF_NEAREST )
  if ( ( t * w->fs - (double) i0 ) < 0.5 ) {
    return w->wf[i0];
  } else return w->wf[i1];

  return -DBL_MAX;
}

// ----------------------------------------------------------------------------

int doublewf_resample( doublewf_t *w2, double fs, doublewf_t *w1, unsigned int mode ) {

  int i = 0;
  if ( ! w1 || ! w2 ) {
    bpm_error( "Invalid pointer arguments in doublewf_resample()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // number of samples
  w2->ns = (int) (w1->ns * fs / w1->fs);
  w2->fs = fs;

  if ( w2->ns > MAX_ALLOWED_NS ) {
    bpm_error( "Maximum allowed number of samples exceeded in doublewf_resample()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  if ( w2->ns < 1 ) {
    bpm_error( "Number of new samples is zero in doublewf_resample()",
	       __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  // fill the new waveform
  for ( i=0; i<w2->ns; i++ ) {
    w2->wf[i] = doublewf_getvalue( w1, (double) i / w2->fs, mode );
  }


  return BPM_SUCCESS;
}
