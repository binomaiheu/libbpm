/**
   @file     bpm_units.h
   @ingroup  general
   @brief    Physical unit definitions for libbpm
*/

#ifndef BPMUNITS_H__
#define BPMUNITS_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <bpm/bpm_defs.h>

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */

#define _cent__ 0.01

#define _Hz__   1.0
#define _kHz__  1.0E+3
#define _MHz__  1.0E+6
#define _GHz__  1.0E+9

#define _sec__  1.0
#define _msec__ 1.0E-3
#define _usec__ 1.0E-6
#define _nsec__ 1.0E-9

#define _eV__   1.0
#define _keV__  1.0E+3
#define _MeV__  1.0E+6
#define _GeV__  1.0E+9

#define  _rad__ 1.0
#define _mrad__ 1.0E-3
#define _urad__ 1.0E-6
#define _nrad__ 1.0E-9

#define _degrees__ (PI/180.)

#define _mC__   1.0E-3
#define _uC__   1.0E-6
#define _nC__   1.0E-9
#define _pC__   1.0E-12

#define _meter__  1.0
#define _mmeter__ 1.0E-3
#define _umeter__ 1.0E-6
#define _nmeter__ 1.0E-9

#define _Volt__  1.0
#define _mVolt__ 1.0E-3
#define _uVolt__ 1.0E-6
#define _nVolt__ 1.0E-9

#define _cLight__ ( 2.99792458E8 * _meter__ / _sec__ )


/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */


/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */


#endif /* #ifndef BPMUNITS_H_LOADED */
/* ================================ end of file ============================= */
