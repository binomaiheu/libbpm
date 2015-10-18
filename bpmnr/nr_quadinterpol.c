/**
   @file
   @ingroup nr
*/

#include "bpm_nr.h"

double nr_quadinterpol( double x, 
			double x1, double x2, double x3,
			double y1, double y2, double y3 ) {

  double d1, d2, d3, d4;
  double a, b, c;

  // if the 3 points are equal, set the interpolated point equal as well
  if( ( fabs( y1 - y2 ) < 1.0e-15 ) &&
      ( fabs( y2 - y3 ) < 1.0e-15 ) ) return y2;
  
  // construct the equation of the parabola from the determinant
  // d1->4 are the 4 3x3 minor determinants
  // d1 * x^2 - d2 * x + d3 * y - d4 = 0
  d1 = x1*y2    + x2*y3    + x3*y1    - x3*y2    - y3*x1    - x2*y1;
  d2 = x1*x1*y2 + x2*x2*y3 + x3*x3*y1 - x3*x3*y2 - y3*x1*x1 - x2*x2*y1;
  d3 = x1*x1*x2 + x2*x2*x3 + x3*x3*x1 - x3*x3*x2 - x3*x1*x1 - x2*x2*x1;
  d4 = x1*x1*x2*y3 + x2*x2*x3*y1 + x1*y2*x3*x3 
     - x3*x3*x2*y1 - x1*x1*x3*y2 - x2*x2*x1*y3;

  if( fabs( d3 ) < 1.e-15 ) {
    //bpm_warning( "Collinearity issue in nr_quadinterpol()", __FILE__, __LINE__ );
    return y2;
  }

  // now calculate a,b,c, don't forget the - signs !
  // y = a * x^2 + b*x + c
  a = - d1 / d3;
  b =   d2 / d3;
  c =   d4 / d3;

  // return interpolated value at the point x
  return a*x*x + b*x + c;
}
