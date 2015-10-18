/** 
   @file      
   @ingroup orbit
   @brief   libbpm orbit generation routines

   This header contains beam oribit generation routines, so this 
   includes also calibration scans etc...
*/

/** @addtogroup orbit */
/** @{ */

#ifndef BPMORBIT_H__
#define BPMORBIT_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <math.h>
#include <bpm/bpm_defs.h>
#include <bpm/bpm_units.h>
#include <bpm/bpm_interface.h>

/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

  /**
     Structure representing a 3-vector, for use in the orbit 
     generation routines
  */
  struct v3 { 
    double x; /**< x-coordinate */
    double y; /**< y-coordinate */
    double z; /**< z-coordinate */
  };
  
  /**
     Structure representing a 3x3-matrix, for use in the orbit 
     generation routines
  */
  struct m33 {
    double e[3][3]; /**< the matrix */
  };



/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */
  
  /**
     Get the bending angle through a rectangular bending magnet
     @param e the particle's charge in units of e, take sign into account !
     @param B the magnetic field in Tesla
     @param l the length of the magnet in meter
     @param p the momentum of the particle in GeV
     @return the bending angle
   */
  EXTERN double get_rbend( double e, double B, double l, double p );


  /**
     Get the bending angle through a sector bending magnet
     @param e the particle's charge in units of e, take sign into account !
     @param B the magnetic field in Tesla
     @param l the sector length of the magnet in meter
     @param p the momentum of the particle in GeV
     @return the bending angle
   */
  EXTERN double get_sbend( double e, double B, double l, double p );
  

  /**
     Get the bunch hit in the local BPM coordinate frame
     @param bunch the bunch structure
     @param bpm   the bpm config
  */
  EXTERN int get_bpmhit( bunchconf_t* bunch, bpmconf_t* bpm );

  
  /**
     Calls get_bpmhit for every bunch in the beam...
     @param beam the beam structure
     @param bpm  the bpm config
  */
  EXTERN int get_bpmhits( beamconf_t* beam, bpmconf_t* bpm );


  /** Copy 3-vector v2 into 3-vector v1 */
  void   v_copy(struct v3 *v1, struct v3 *v2);

  /** Return the magnitude of 3-vector v1 */
  double v_mag(struct v3 *v1);

  /** Scale 3-vector v1 with factor dscale */
  void   v_scale(struct v3 *v1, double dscale);

  /** Normalise 3-vector v1 to unit vector */
  void   v_norm(struct v3 *v1);

  /** Multiply matrix m1 with 3-vector v1 : m1.v1, result is in v1 */
  void   v_matmult(struct m33 *m1, struct v3 *v1);

  /** Add two 3-vectors v1 and v2, result is in v1 */
  void   v_add(struct v3 *v1, struct v3 *v2);

  /** Subtract 3-vectors v1 - v2, result is in v1*/
  void   v_sub(struct v3 *v1, struct v3 *v2);

  /** Return Scalar product of 3-vectors v1 and v2 */
  double v_dot(struct v3 *v1, struct v3 *v2);

  /** Return the vector product of 3 vectors v1 x v2, result is in v1 */
  void   v_cross(struct v3 *v1, struct v3 *v2);

  /** Print the 3-vector to stdout */
  void   v_print(struct v3 *v1);

  /** Create rotation 3x3 matrix with the 3 euler angles alpha, beta and gamma, result in m1 */
  void   m_rotmat(struct m33 *m1, double alpha, double beta, double gamma);

  /** 3x3 Matrix multiplication m1.m2, result in m*/
  void   m_matmult(struct m33 *m, struct m33 *m1, struct m33 *m2);

  /** 3x3 Matrix addition m1+m2, result in m1 */
  void   m_matadd(struct m33 *m1, struct m33 *m2);

  /** Print 3x3 matrix m1 to stdout */
  void   m_print(struct m33 *m1);
  
#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMORBIT_H__ */
/* @} */
/* ================================ end of file ============================= */
