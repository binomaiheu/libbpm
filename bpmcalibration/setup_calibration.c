/**
   @file
   @ingroup calib
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_calibration.h>

int setup_calibration( bpmconf_t *cnf, bpmproc_t *proc, int npulses, int startpulse, 
		       int stoppulse, double angle, double startpos, double endpos, 
		       int num_steps, bunchconf_t *bunch ) {

  if ( ! proc || ! bunch || ! cnf ) {
    bpm_error( "Invalid pointer arguments in setup_calibration(...)",
               __FILE__, __LINE__ );
    return BPM_FAILURE;
  }

  int ns = 5, nerr = 0, nstep, i;
  double mean = 0, stddev = 0, val, prev_mean, prev_stddev, cnt = 0;
  double curr_beam_pos = startpos, stepsize = (endpos - startpos) / (double) num_steps;

  for ( nstep = 0; nstep < num_steps; nstep++) { 
    
    while ( nerr < 3 ){
      
      mean = 0;
      stddev = 0;
      cnt = 0;
      
      // Calculate mean and std. dev. from first ns pulses
      for (i = startpulse; i < ns + startpulse; i++ ) {
	val = sqrt( pow( ( proc[i].ddc_I - proc[stoppulse].ddc_I ) , 2) + 
		    pow( ( proc[i].ddc_Q - proc[stoppulse].ddc_Q ) , 2) );
	
	if ( fabs( val - prev_mean ) < ( prev_stddev * 3.) ) {
	  mean += sqrt( pow( ( proc[i].ddc_I - proc[stoppulse].ddc_I ) , 2) + 
			pow( ( proc[i].ddc_Q - proc[stoppulse].ddc_Q ) , 2) );
	}
	
	cnt++;
      }
      
      mean /= cnt;
      
      for (i = startpulse; i < ns + startpulse; i++ ) {
	val = sqrt( pow( ( proc[i].ddc_I - proc[stoppulse].ddc_I ) , 2) + 
		    pow( ( proc[i].ddc_Q - proc[stoppulse].ddc_Q ) , 2) );
	
	if ( fabs( val - prev_mean ) < ( prev_stddev * 3.) ) {
	  stddev += (val - mean) * (val - mean);
	}
      }
      
      stddev = sqrt( stddev / cnt );
      
      prev_mean = mean;
      prev_stddev = stddev;
      
      // Add the next event?
      val = sqrt( pow( ( proc[ns + startpulse].ddc_I - proc[stoppulse].ddc_I ) , 2) + 
		  pow( ( proc[ns + startpulse].ddc_Q - proc[stoppulse].ddc_Q ) , 2) );
      
      if ( fabs( val - mean ) > ( stddev * 3.) ) nerr++;
      
      ns++;
    }
    
    // Add the beam position to the beamconf array
    for (i = startpulse; i < ns + startpulse; i++ ) {
      
      bunch[i].bpmposition[0] = cos( angle ) * curr_beam_pos;
      bunch[i].bpmposition[1] = sin( angle ) * curr_beam_pos;
    }
    
    curr_beam_pos += stepsize;
    startpulse += ns; 
  }

  return BPM_SUCCESS;
}


