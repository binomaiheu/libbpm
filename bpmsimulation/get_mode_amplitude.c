/**
   @file
   @ingroup sim
*/

#include <bpm/bpm_simulation.h>
#include <math.h>

complex_t get_mode_amplitude( bpmconf_t *bpm, bpmmode_t *mode, bunchconf_t *bunch ) {

  double k, L, sigma;
  complex_t z;
  int idx = ( mode->polarisation == horiz ? 0 : 1);

  if ( mode->order == 0 ) {
    // monopole mode, sensititvity in Volt/C
    z.re = mode->sensitivity * bunch->charge;
    z.im = 0.;

  } else if ( mode->order == 1 ) {
    // dipole mode, sensitivity in Volt/m/C
    // need to include the bunch length effects properly!!!
    z.re = mode->sensitivity * bunch->bpmposition[idx] * bunch->charge;

    // calulcate :
    // beam->bunchlength, bpm->cav_length, beam->bpmslope[idx] + beam->bpmtilt[idx]
    k = 2*PI / _cLight__ * mode->frequency;
    L = bpm->cav_length;
    sigma = bunch->length;

    if ( L != 0 && bunch->bpmslope[mode->polarisation] !=0 ) {

      z.im = mode->sensitivity * bunch->charge * 
             ( 1/k * ( 1 - (k*L/2) / sin (k*L/2) ) *
             bunch->bpmslope[mode->polarisation]
	     + 2 * sigma * sin (k*sigma/2) ); 

    } else {

      z.im = 0;

    }

  } else if (mode->order == 2 ) {
    // quadrupole mode
    z.re = 0.;
    z.im = 0.;
    bpm_warning( "Quadrupole modes are not implemented yet in libbpm...", __FILE__, __LINE__ );

  } else {
    bpm_warning( "HOM (O(>2) modes are not implemented yet in libbpm...", __FILE__, __LINE__ );

  }

  return z;
}
