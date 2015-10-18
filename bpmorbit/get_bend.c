/**
   get_rbend.c
*/

#include <bpm/bpm_orbit.h>

// B: Magnetic field in Tesla
// l: Length in meter
// p: Momentum 
// e: charge in units of e, include sign !!

double get_rbend( double e, double B, double l, double p ) {
  // l is rectangular magnet length
  return asin( e * _cLight__ * B * l /  p );
}

double get_sbend( double e, double B, double l, double p ) {
  // l is sector length
  return e * _cLight__ * B * l /  p;
}
