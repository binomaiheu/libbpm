/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

double gsl_blas_dnrm2 (const gsl_vector * X)
{
    /**
    * 
    * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
  
  return cblas_dnrm2 ((int)(X->size), X->data, (int)(X->stride));
}

double cblas_dnrm2(const int N, const double *X, const int incX) {


  double scale = 0.0;
  double ssq = 1.0;
  int i;
  int ix = 0;

  if (N <= 0 || incX <= 0) {
    return 0;
  } else if (N == 1) {
    return fabs(X[0]);
  }

  for (i = 0; i < N; i++) {
    const double x = X[ix];

    if (x != 0.0) {
      const double ax = fabs(x);

      if (scale < ax) {
        ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
        scale = ax;
      } else {
        ssq += (ax / scale) * (ax / scale);
      }
    }

    ix += incX;
  }

  return scale * sqrt(ssq);
}

void gsl_blas_dscal (double alpha, gsl_vector * X)
{
  cblas_dscal ((int) (X->size), alpha, X->data, (int) (X->stride));
}

void cblas_dscal (const int N, const double alpha, double *X, const int incX)
{
  int i;
  int ix;

  if (incX <= 0) {
    return;
  }

  ix = OFFSET(N, incX);

  for (i = 0; i < N; i++) {
    X[ix] *= alpha;
    ix += incX;
  }
}

int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A,
		    const gsl_vector * X, double beta, gsl_vector * Y)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if ((TransA == CblasNoTrans && N == X->size && M == Y->size)
      || (TransA == CblasTrans && M == X->size && N == Y->size))
    {
      cblas_dgemv (CblasRowMajor, TransA, (int) M, (int) N, alpha, A->data,
                   (int) (A->tda), X->data, (int) (X->stride), beta, Y->data,
                   (int) (Y->stride));
      return BPM_SUCCESS;
    }
  else
    {
      bpm_error("invalid length in gsl_blas_dgemv(..)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
}

void
cblas_dgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
             const int M, const int N, const double alpha, const double *A,
             const int lda, const double *X, const int incX,
             const double beta, double *Y, const int incY)
{
  int i, j;
  int lenX, lenY;

  const int Trans = (TransA != CblasConjTrans) ? TransA : CblasTrans;

  if (M == 0 || N == 0)
    return;

  if (alpha == 0.0 && beta == 1.0)
    return;

  if (Trans == CblasNoTrans) {
    lenX = N;
    lenY = M;
  } else {
    lenX = M;
    lenY = N;
  }

  /* form  y := beta*y */
  if (beta == 0.0) {
    int iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] = 0.0;
      iy += incY;
    }
  } else if (beta != 1.0) {
    int iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      Y[iy] *= beta;
      iy += incY;
    }
  }

  if (alpha == 0.0)
    return;

  if ((order == CblasRowMajor && Trans == CblasNoTrans)
      || (order == CblasColMajor && Trans == CblasTrans)) {
    /* form  y := alpha*A*x + y */
    int iy = OFFSET(lenY, incY);
    for (i = 0; i < lenY; i++) {
      double temp = 0.0;
      int ix = OFFSET(lenX, incX);
      for (j = 0; j < lenX; j++) {
        temp += X[ix] * A[lda * i + j];
        ix += incX;
      }
      Y[iy] += alpha * temp;
      iy += incY;
    }
  } else if ((order == CblasRowMajor && Trans == CblasTrans)
             || (order == CblasColMajor && Trans == CblasNoTrans)) {
    /* form  y := alpha*A'*x + y */
    int ix = OFFSET(lenX, incX);
    for (j = 0; j < lenX; j++) {
      const double temp = alpha * X[ix];
      if (temp != 0.0) {
        int iy = OFFSET(lenY, incY);
        for (i = 0; i < lenY; i++) {
          Y[iy] += temp * A[lda * j + i];
          iy += incY;
        }
      }
      ix += incX;
    }
  } else {
    bpm_error("unrecognised operation in cblas_dgemv(..)",
	      __FILE__, __LINE__);
    return;
  }
}
