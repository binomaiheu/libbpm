/** 
    @file
    @ingroup nr
    @brief   libbpm numerical helper routines

    Header file containing the numerical recipies and GNU Scientific
    Library routines used in the library.
*/

/** @addtogroup nr */
/** @{ */

#ifndef BPMNR_H__
#define BPMNR_H__

/* -----------------------------------------------------------------------------
   includes
   -------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include <bpm/bpm_defs.h>

/* -----------------------------------------------------------------------------
   macro definitions
   -------------------------------------------------------------------------- */
#define GCF_ITMAX  100
#define GCF_FPMIN  1.0e-30
#define GCF_EPS    3.0e-7

#define GSER_EPS   3.0e-7
#define GSER_ITMAX 100

#define RAN1_IA 16807
#define RAN1_IM 2147483647
#define RAN1_AM (1.0/RAN1_IM)
#define RAN1_IQ 127773
#define RAN1_IR 2836
#define RAN1_NTAB 32
#define RAN1_NDIV (1+(RAN1_IM-1)/RAN1_NTAB)
#define RAN1_EPS 1.2e-7
#define RAN1_RNMX (1.0-RAN1_EPS)

/**
 Block size for cache-friendly matrix-matrix multiply. It should be
 such that __BLOCKSZ__^2*sizeof(LM_REAL) is smaller than the CPU (L1)
 data cache size. Notice that a value of 32 when LM_REAL=double 
 assumes an 8Kb L1 data cache (32*32*8=8K). This is a concervative 
 choice since newer Pentium 4s have a L1 data cache of size 16K, 
 capable of holding up to 45x45 double blocks.
*/
#define __LM_BLOCKSZ__ 32
#define __LM_BLOCKSZ__SQ    (__LM_BLOCKSZ__)*(__LM_BLOCKSZ__)

#define LINSOLVERS_RETAIN_MEMORY
#ifdef LINSOLVERS_RETAIN_MEMORY
#define __LM_STATIC__ static
#else
#define __LM_STATIC__ /* empty */
#endif /* LINSOLVERS_RETAIN_MEMORY */

#define FABS(x) (((x)>=0.0)? (x) : -(x))
#define CNST(x) (x)
#define _LM_POW_ CNST(2.1)

/**
   Work array size for LM with & without jacobian, should be multiplied by sizeof(double)
   or sizeof(float) to be converted to bytes
*/
#define LM_DER_WORKSZ(npar, nmeas) (2*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))
/** see LM_DER_WORKSZ */
#define LM_DIF_WORKSZ(npar, nmeas) (3*(nmeas) + 4*(npar) + (nmeas)*(npar) + (npar)*(npar))

#define LM_EPSILON       1E-12
#define LM_ONE_THIRD     0.3333333334 /* 1.0/3.0 */

#define LM_OPTS_SZ    	 5 /* max(4, 5) */
#define LM_INFO_SZ    	 9
#define LM_INIT_MU    	 1E-03
#define LM_STOP_THRESH	 1E-17
#define LM_DIFF_DELTA    1E-06

#define NR_FFTFORWARD     1  /**< Perform forward FFT in nr_four */
#define NR_FFTBACKWARD   -1  /**< Perform backward FFT in nr_four */

/** find the median of 3 numbers */
#define __LM_MEDIAN3(a, b, c) ( ((a) >= (b))?\
        ( ((c) >= (a))? (a) : ( ((c) <= (b))? (b) : (c) ) ) : \
        ( ((c) >= (b))? (b) : ( ((c) <= (a))? (a) : (c) ) ) )

/* Some GSL defines */
#define NULL_VECTOR {0, 0, 0, 0, 0}
#define NULL_VECTOR_VIEW {{0, 0, 0, 0, 0}}
#define NULL_MATRIX {0, 0, 0, 0, 0, 0}
#define NULL_MATRIX_VIEW {{0, 0, 0, 0, 0, 0}}
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define OFFSET(N, incX) ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))


/* -----------------------------------------------------------------------------
   typedefs, struct, enums and other declarations
   -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif

  enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
  enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
  typedef  enum CBLAS_TRANSPOSE   CBLAS_TRANSPOSE_t;
  
  /**
     structure needed for levenberg marquard minimisation 
  */
  struct lm_fstate{ 
    int n, *nfev;
    double *hx, *x;
    void *adata;
  };

  /*
    Here are some GSL structures and the licence info. Please don't sue us!
    * 
    * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman, Brian Gough
    * 
    * This program is free software; you can redistribute it and/or modify
    * it under the terms of the GNU General Public License as published by
    * the Free Software Foundation; either version 2 of the License, or (at
    * your option) any later version.
    * 
    * This program is distributed in the hope that it will be useful, but
    * WITHOUT ANY WARRANTY; without even the implied warranty of
    * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    * General Public License for more details.
    * 
    * You should have received a copy of the GNU General Public License
    * along with this program; if not, write to the Free Software
    * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
    */
  
/* --------------------------------------
   Blocks */
  struct gsl_block_struct
  {
    size_t size;
    double *data;
  };
  
  typedef struct gsl_block_struct gsl_block;
  
/* --------------------------------------
   Matrices */
  typedef struct 
  {
    size_t size1;
    size_t size2;
    size_t tda;
    double * data;
    gsl_block * block;
    int owner;
  } gsl_matrix;
  
  typedef struct
  {
    gsl_matrix matrix;
  } 
    _gsl_matrix_view;
  
  typedef _gsl_matrix_view gsl_matrix_view;
  
/* --------------------------------------
   Vectors */
  typedef struct 
  {
    size_t size;
    size_t stride;
    double *data;
    gsl_block *block;
    int owner;
  } 
    gsl_vector;
  
  typedef struct
  {
    gsl_vector vector;
  }
    _gsl_vector_view;
  
  typedef _gsl_vector_view gsl_vector_view;
  
  typedef struct
  {
    gsl_vector vector;
  } 
    _gsl_vector_const_view;
  
  typedef const _gsl_vector_const_view gsl_vector_const_view;
  
  
  /**
     Structure and typedef for complex numbers used in the bpmdsp module
  */
  typedef struct {
    double re;
    double im;
  } complex_t;

/* -----------------------------------------------------------------------------
   function prototypes and declarations
   -------------------------------------------------------------------------- */
  
  EXTERN double nr_gammln( double xx );
  EXTERN double nr_gammq( double a, double x );
  EXTERN int nr_gcf( double *gammcf, double a, double x, double *gln );
  EXTERN int nr_gser( double *gamser, double a, double x, double *gln );
  EXTERN int nr_fit ( double *x, double y[], int ndata, double sig[], 
		      int mwt, double *a, double *b, 
		      double *siga, double *sigb, double *chi2, double *q );
  EXTERN int nr_is_pow2( unsigned long n );
  EXTERN int nr_four1( double data[], unsigned long nn, int isign );
  EXTERN int nr_realft( double data[], unsigned long n, int isign );
  EXTERN double nr_ran1( long *idum );
  EXTERN int nr_seed( long seed );
  EXTERN double nr_ranuniform( double lower, double upper );
  EXTERN double nr_rangauss( double mean, double std_dev );
  EXTERN long bpm_rseed;

  /* the levenberg marquardt routines from nr_levmar.c */
  /* unconstrained minimization with known jacobian */
  EXTERN int nr_lmder( void (*func)(double *p, double *hx, int m, int n, void *adata),
		       void (*jacf)(double *p, double *j, int m, int n, void *adata),
		       double *p, double *x, int m, int n, int itmax, double *opts,
		       double *info, double *work, double *covar, void *adata );
  
  /* unconstrained minimization with unknown jacobian */
  EXTERN int nr_lmdif( void (*func)(double *p, double *hx, int m, int n, void *adata),
		       double *p, double *x, int m, int n, int itmax, double *opts,
		       double *info, double *work, double *covar, void *adata );

  /* box-constrained minimization with known jacobian */
  EXTERN int nr_lmder_bc( void (*func)(double *p, double *hx, int m, int n, void *adata),
			  void (*jacf)(double *p, double *j, int m, int n, void *adata),  
			  double *p, double *x, int m, int n, double *lb, double *ub, int itmax,
			  double *opts, double *info, double *work, double *covar, void *adata );

  /* box-constrained minimization with unknown jacobian */
  EXTERN int nr_lmdif_bc( void (*func)(double *p, double *hx, int m, int n, void *adata),
			  double *p, double *x, int m, int n, double *lb, double *ub, int itmax, 
			  double *opts, double *info, double *work, double *covar, void *adata );

  /* jacobian verification */
  EXTERN void nr_lmchkjac( void (*func)(double *p, double *hx, int m, int n, void *adata),
			   void (*jacf)(double *p, double *j, int m, int n, void *adata),
			   double *p, int m, int n, void *adata, double *err );
  
  /* covariance matrix */
  EXTERN int nr_lmcovar(double *JtJ, double *C, double sumsq, int m, int n );

  /* lm helper routine, LU solver for matrix inversion */
  EXTERN int nr_ax_eq_b_LU(double *A, double *B, double *x, int n);
  EXTERN void nr_trans_mat_mat_mult(double *a, double *b, int n, int m);
  EXTERN void nr_fdif_forw_jac_approx( void (*func)(double *p, double *hx, 
						    int m, int n, void *adata),
				       double *p, double *hx, double *hxx, double delta,
				       double *jac, int m, int n, void *adata );
  
  EXTERN void nr_fdif_cent_jac_approx(void (*func)(double *p, double *hx, 
						   int m, int n, void *adata),
				      double *p, double *hxm, double *hxp, double delta,
				      double *jac, int m, int n, void *adata );

  EXTERN double nr_median( int n, double *arr );

  EXTERN double nr_select( int k, int n, double *org_arr );
  
  /* 
     Here are a load of GNU Scientific Library functions. See above for the license info.
  */

  /* ----------------------------------------
     Matrix related */
  EXTERN gsl_matrix * gsl_matrix_calloc (const size_t n1, const size_t n2);
  EXTERN _gsl_vector_view gsl_matrix_column (gsl_matrix * m, const size_t i);
  EXTERN _gsl_matrix_view gsl_matrix_submatrix (gsl_matrix * m, 
                            const size_t i, const size_t j, 
                            const size_t n1, const size_t n2);
  EXTERN double gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j); 
  EXTERN void gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x);
  EXTERN int gsl_matrix_swap_columns(gsl_matrix * m, const size_t i, const size_t j);
  EXTERN gsl_matrix * gsl_matrix_alloc (const size_t n1, const size_t n2);
  EXTERN _gsl_vector_const_view gsl_matrix_const_row (const gsl_matrix * m, const size_t i);
  EXTERN _gsl_vector_view gsl_matrix_row (gsl_matrix * m, const size_t i);
  EXTERN _gsl_vector_const_view gsl_matrix_const_column (const gsl_matrix * m, const size_t j);
  EXTERN void gsl_matrix_set_identity (gsl_matrix * m);
  
  /* ----------------------------------------
     Vector related */
  EXTERN gsl_vector *gsl_vector_calloc (const size_t n);
  EXTERN _gsl_vector_view gsl_vector_subvector (gsl_vector *v, size_t offset, size_t n);
  EXTERN double gsl_vector_get (const gsl_vector * v, const size_t i);
  EXTERN void gsl_vector_set (gsl_vector * v, const size_t i, double x);
  EXTERN int gsl_vector_swap_elements (gsl_vector * v, const size_t i, const size_t j);
  EXTERN _gsl_vector_const_view gsl_vector_const_subvector (const gsl_vector *v, size_t i, 
							    size_t n);
  EXTERN void gsl_vector_free (gsl_vector * v);

  /* ----------------------------------------
     linalg related */
  EXTERN int gsl_linalg_SV_solve (const gsl_matrix * U, const gsl_matrix * Q,
				  const gsl_vector * S, const gsl_vector * b,
				  gsl_vector * x);
  EXTERN int gsl_linalg_bidiag_unpack (const gsl_matrix * A, const gsl_vector * tau_U, 
				       gsl_matrix * U, const gsl_vector * tau_V,
				       gsl_matrix * V, gsl_vector * diag, 
				       gsl_vector * superdiag);
  EXTERN int gsl_linalg_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A);
  EXTERN int gsl_linalg_bidiag_unpack2 (gsl_matrix * A, gsl_vector * tau_U, 
					gsl_vector * tau_V, gsl_matrix * V);
  EXTERN int gsl_linalg_householder_hm1 (double tau, gsl_matrix * A);
  EXTERN void create_givens (const double a, const double b, double *c, double *s);
  EXTERN double gsl_linalg_householder_transform (gsl_vector * v);
  EXTERN int gsl_linalg_householder_mh (double tau, const gsl_vector * v, gsl_matrix * A);
  EXTERN int gsl_linalg_SV_solve (const gsl_matrix * U,
				  const gsl_matrix * V,
				  const gsl_vector * S,
				  const gsl_vector * b, gsl_vector * x);
  EXTERN void chop_small_elements (gsl_vector * d, gsl_vector * f);
  EXTERN void qrstep (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V);
  EXTERN double trailing_eigenvalue (const gsl_vector * d, const gsl_vector * f);
  EXTERN void create_schur (double d0, double f0, double d1, double * c, double * s);
  EXTERN void svd2 (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V);
  EXTERN void chase_out_intermediate_zero (gsl_vector * d, gsl_vector * f, gsl_matrix * U, 
					   size_t k0);
  EXTERN void chase_out_trailing_zero (gsl_vector * d, gsl_vector * f, gsl_matrix * V);
  EXTERN int gsl_isnan (const double x);

  /* ----------------------------------------
     blas/cblas related */
  EXTERN double gsl_blas_dnrm2 (const gsl_vector * X);
  EXTERN double cblas_dnrm2(const int N, const double *X, const int incX);
  EXTERN void gsl_blas_dscal  (double alpha, gsl_vector * X);
  EXTERN void cblas_dscal(const int N, const double alpha, double *X, const int incX);
  EXTERN void cblas_dgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
			   const int M, const int N, const double alpha, const double *A,
			   const int lda, const double *X, const int incX,
			   const double beta, double *Y, const int incY);

  /* ----------------------------------------
     block related */
  EXTERN gsl_block *gsl_block_alloc (const size_t n);
  EXTERN void gsl_block_free (gsl_block * b);



  /* Some complex functions */
  EXTERN complex_t complex( double re, double im );
  EXTERN double    c_real( complex_t z );
  EXTERN double    c_imag( complex_t z );
  EXTERN complex_t c_conj( complex_t z );
  EXTERN complex_t c_neg( complex_t z );
  EXTERN complex_t c_sum( complex_t z1, complex_t z2 );
  EXTERN complex_t c_diff( complex_t z1, complex_t z2 );
  EXTERN complex_t c_mult( complex_t z1, complex_t z2 );
  EXTERN complex_t c_div( complex_t z1, complex_t z2 );
  EXTERN complex_t c_scale( double r, complex_t z );
  EXTERN complex_t c_sqr( complex_t z );
  EXTERN complex_t c_sqrt( complex_t z );
  EXTERN double    c_norm2( complex_t z );
  EXTERN double    c_abs( complex_t z );
  EXTERN double    c_abs2( complex_t z );
  EXTERN double    c_arg( complex_t z );
  EXTERN complex_t c_exp( complex_t z );
  EXTERN int       c_isequal( complex_t z1, complex_t z2 );


  /** Parabolic (quadratic) interpolation routine, give 3 points (x1,y1),
      (x2,y2) and (x3,y3) and a value x which needs to be interpolated. 
      The function returns y, which is the value of a parabola at point x
      defined by the 3 points given */
  EXTERN double nr_quadinterpol( double x, 
				 double x1, double x2, double x3,
				 double y1, double y2, double y3 );
  
  
  /** The normalised sinc(x) function */
  EXTERN double sinc( double x );

  /** The Lanczos kernel */
  EXTERN double lanczos( double x, int a );

  /** Rounds a value to nearest integers, voids the need for -std=c99 in the compilation */
  EXTERN double dround( double x );


#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMNR_H__ */
/*@}*/
/* ================================ end of file ============================= */
