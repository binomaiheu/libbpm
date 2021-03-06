/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

int gsl_linalg_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /**
    applies a householder transformation v,tau to matrix m 

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

  if (tau == 0.0)
    {
      return BPM_SUCCESS;
    }

#ifdef USE_BLAS
  {
    gsl_vector_const_view v1 = gsl_vector_const_subvector (v, 1, v->size - 1);
    gsl_matrix_view A1 = gsl_matrix_submatrix (A, 1, 0, A->size1 - 1, A->size2);
    size_t j;

    for (j = 0; j < A->size2; j++)
      {
        double wj = 0.0;
        gsl_vector_view A1j = gsl_matrix_column(&A1.matrix, j);
        gsl_blas_ddot (&A1j.vector, &v1.vector, &wj);
        wj += gsl_matrix_get(A,0,j);

        {
          double A0j = gsl_matrix_get (A, 0, j);
          gsl_matrix_set (A, 0, j, A0j - tau *  wj);
        }

        gsl_blas_daxpy (-tau * wj, &v1.vector, &A1j.vector);
      }
  }
#else
  {
    size_t i, j;
    
    for (j = 0; j < A->size2; j++)
      {
        /* Compute wj = Akj vk */
        
        double wj = gsl_matrix_get(A,0,j);  
        
        for (i = 1; i < A->size1; i++)  /* note, computed for v(0) = 1 above */
          {
            wj += gsl_matrix_get(A,i,j) * gsl_vector_get(v,i);
          }
        
        /* Aij = Aij - tau vi wj */
        
        /* i = 0 */
        {
          double A0j = gsl_matrix_get (A, 0, j);
          gsl_matrix_set (A, 0, j, A0j - tau *  wj);
        }
        
        /* i = 1 .. M-1 */
        
        for (i = 1; i < A->size1; i++)
          {
            double Aij = gsl_matrix_get (A, i, j);
            double vi = gsl_vector_get (v, i);
            gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
          }
      }
  }
#endif
    
  return BPM_SUCCESS;
}

int gsl_linalg_householder_hm1 (double tau, gsl_matrix * A)
{
  /**
    applies a householder transformation v,tau to a matrix being
    build up from the identity matrix, using the first column of A as
    a householder vector
  */

  if (tau == 0)
    {
      size_t i,j;

      gsl_matrix_set (A, 0, 0, 1.0);
      
      for (j = 1; j < A->size2; j++)
        {
          gsl_matrix_set (A, 0, j, 0.0);
        }

      for (i = 1; i < A->size1; i++)
        {
          gsl_matrix_set (A, i, 0, 0.0);
        }

      return BPM_SUCCESS;
    }

  /* w = A' v */

#ifdef USE_BLAS
  {
    gsl_matrix_view A1 = gsl_matrix_submatrix (A, 1, 0, A->size1 - 1, A->size2);
    gsl_vector_view v1 = gsl_matrix_column (&A1.matrix, 0);
    size_t j;

    for (j = 1; j < A->size2; j++)
      {
        double wj = 0.0;   /* A0j * v0 */
        
        gsl_vector_view A1j = gsl_matrix_column(&A1.matrix, j);
        gsl_blas_ddot (&A1j.vector, &v1.vector, &wj);

        /* A = A - tau v w' */
        
        gsl_matrix_set (A, 0, j, - tau *  wj);
        
        gsl_blas_daxpy(-tau*wj, &v1.vector, &A1j.vector);
      }

    gsl_blas_dscal(-tau, &v1.vector);
    
    gsl_matrix_set (A, 0, 0, 1.0 - tau);
  }
#else
  {
    size_t i, j;
    
    for (j = 1; j < A->size2; j++)
      {
        double wj = 0.0;   /* A0j * v0 */
        
        for (i = 1; i < A->size1; i++)
          {
            double vi = gsl_matrix_get(A, i, 0);
            wj += gsl_matrix_get(A,i,j) * vi;
          }
        
        /* A = A - tau v w' */
        
        gsl_matrix_set (A, 0, j, - tau *  wj);
        
        for (i = 1; i < A->size1; i++)
          {
            double vi = gsl_matrix_get (A, i, 0);
            double Aij = gsl_matrix_get (A, i, j);
            gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
          }
      }
    
    for (i = 1; i < A->size1; i++)
      {
        double vi = gsl_matrix_get(A, i, 0);
        gsl_matrix_set(A, i, 0, -tau * vi);
      }
    
    gsl_matrix_set (A, 0, 0, 1.0 - tau);
  }
#endif

  return BPM_SUCCESS;
}

void create_givens (const double a, const double b, double *c, double *s)
{
  if (b == 0)
    {
      *c = 1;
      *s = 0;
    }
  else if (fabs (b) > fabs (a))
    {
      double t = -a / b;
      double s1 = 1.0 / sqrt (1 + t * t);
      *s = s1;
      *c = s1 * t;
    }
  else
    {
      double t = -b / a;
      double c1 = 1.0 / sqrt (1 + t * t);
      *c = c1;
      *s = c1 * t;
    }
}


int gsl_linalg_bidiag_decomp (gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V)  
{
  if (A->size1 < A->size2)
    {
      bpm_error("bidiagonal decomposition requires M>=N in gsl_linalg_bidag_decomp(...)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else if (tau_U->size  != A->size2)
    {
      bpm_error("size of tau_U must be N in gsl_linalg_bidag_decomp(...)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else if (tau_V->size + 1 != A->size2)
    {
      bpm_error("size of tau_V must be (N - 1) in gsl_linalg_bidag_decomp(...)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else
    {
      const size_t M = A->size1;
      const size_t N = A->size2;
      size_t i;
  
      for (i = 0 ; i < N; i++)
        {
          /* Apply Householder transformation to current column */
          
          {
            gsl_vector_view c = gsl_matrix_column (A, i);
            gsl_vector_view v = gsl_vector_subvector (&c.vector, i, M - i);
            double tau_i = gsl_linalg_householder_transform (&v.vector);
            
            /* Apply the transformation to the remaining columns */
            
            if (i + 1 < N)
              {
                gsl_matrix_view m = 
                  gsl_matrix_submatrix (A, i, i + 1, M - i, N - (i + 1));
                gsl_linalg_householder_hm (tau_i, &v.vector, &m.matrix);
              }

            gsl_vector_set (tau_U, i, tau_i);            

          }

          /* Apply Householder transformation to current row */
          
          if (i + 1 < N)
            {
              gsl_vector_view r = gsl_matrix_row (A, i);
              gsl_vector_view v = gsl_vector_subvector (&r.vector, i + 1, N - (i + 1));
              double tau_i = gsl_linalg_householder_transform (&v.vector);
              
              /* Apply the transformation to the remaining rows */
              
              if (i + 1 < M)
                {
                  gsl_matrix_view m = 
                    gsl_matrix_submatrix (A, i+1, i+1, M - (i+1), N - (i+1));
                  gsl_linalg_householder_mh (tau_i, &v.vector, &m.matrix);
                }

              gsl_vector_set (tau_V, i, tau_i);
            }
        }
    }
        
  return BPM_SUCCESS;
}

double gsl_linalg_householder_transform (gsl_vector * v)
{
  /**
    replace v[0:n-1] with a householder vector (v[0:n-1]) and
    coefficient tau that annihilate v[1:n-1]
  */

  const size_t n = v->size ;

  if (n == 1)
    {
      return 0.0; /* tau = 0 */
    }
  else
    { 
      double alpha, beta, tau ;
      
      gsl_vector_view x = gsl_vector_subvector (v, 1, n - 1) ; 
      
      double xnorm = gsl_blas_dnrm2 (&x.vector);
      
      if (xnorm == 0) 
        {
          return 0.0; /* tau = 0 */
        }
      
      alpha = gsl_vector_get (v, 0) ;
      beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
      tau = (beta - alpha) / beta ;
      
      gsl_blas_dscal (1.0 / (alpha - beta), &x.vector);
      gsl_vector_set (v, 0, beta) ;
      
      return tau;
    }
}

int gsl_linalg_householder_mh (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /**
    applies a householder transformation v,tau to matrix m from the
    right hand side in order to zero out rows
  */

  if (tau == 0)
    return BPM_SUCCESS;

  /* A = A - tau w v' */

#ifdef USE_BLAS
  {
    gsl_vector_const_view v1 = gsl_vector_const_subvector (v, 1, v->size - 1);
    gsl_matrix_view A1 = gsl_matrix_submatrix (A, 0, 1, A->size1, A->size2-1);
    size_t i;

    for (i = 0; i < A->size1; i++)
      {
        double wi = 0.0;
        gsl_vector_view A1i = gsl_matrix_row(&A1.matrix, i);
        gsl_blas_ddot (&A1i.vector, &v1.vector, &wi);
        wi += gsl_matrix_get(A,i,0);  
        
        {
          double Ai0 = gsl_matrix_get (A, i, 0);
          gsl_matrix_set (A, i, 0, Ai0 - tau *  wi);
        }
        
        gsl_blas_daxpy(-tau * wi, &v1.vector, &A1i.vector);
      }
  }
#else
  {
    size_t i, j;
    
    for (i = 0; i < A->size1; i++)
      {
        double wi = gsl_matrix_get(A,i,0);  
        
        for (j = 1; j < A->size2; j++)  /* note, computed for v(0) = 1 above */
          {
            wi += gsl_matrix_get(A,i,j) * gsl_vector_get(v,j);
          }
        
        /* j = 0 */
        
        {
          double Ai0 = gsl_matrix_get (A, i, 0);
          gsl_matrix_set (A, i, 0, Ai0 - tau *  wi);
        }
        
        /* j = 1 .. N-1 */
        
        for (j = 1; j < A->size2; j++) 
          {
            double vj = gsl_vector_get (v, j);
            double Aij = gsl_matrix_get (A, i, j);
            gsl_matrix_set (A, i, j, Aij - tau * wi * vj);
          }
      }
  }
#endif
    
  return BPM_SUCCESS;
}

int gsl_linalg_SV_solve (const gsl_matrix * U,
			 const gsl_matrix * V,
			 const gsl_vector * S,
			 const gsl_vector * b, gsl_vector * x)
{
  if (U->size1 != b->size)
    {
      bpm_error("first dimension of matrix U must size of vector b in gsl_linalg_SV_solve(..)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else if (U->size2 != S->size)
    {
      bpm_error("length of vector S must match second dimension of matrix U in gsl_linalg_SV_solve(..)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else if (V->size1 != V->size2)
    {
      bpm_error("matrix V must be square in gsl_linalg_SV_solve(..)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else if (S->size != V->size1)
    {
      bpm_error("length of vector S must match size of matrix V in gsl_linalg_SV_solve(..)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else if (V->size2 != x->size)
    {
      bpm_error("size of matrix V must match size of vector x in gsl_linalg_SV_solve(..)",
		__FILE__, __LINE__);
      return BPM_SUCCESS;
    }
  else
    {
      const size_t N = U->size2;
      size_t i;

      gsl_vector *w = gsl_vector_calloc (N);

      gsl_blas_dgemv (CblasTrans, 1.0, U, b, 0.0, w);

      for (i = 0; i < N; i++)
        {
          double wi = gsl_vector_get (w, i);
          double alpha = gsl_vector_get (S, i);
          if (alpha != 0)
            alpha = 1.0 / alpha;
          gsl_vector_set (w, i, alpha * wi);
        }

      gsl_blas_dgemv (CblasNoTrans, 1.0, V, w, 0.0, x);

      gsl_vector_free (w);

      return BPM_SUCCESS;
    }
}

int gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}

void chop_small_elements (gsl_vector * d, gsl_vector * f)
{
  const size_t N = d->size;
  double d_i = gsl_vector_get (d, 0);

  size_t i;

  for (i = 0; i < N - 1; i++)
    {
      double f_i = gsl_vector_get (f, i);
      double d_ip1 = gsl_vector_get (d, i + 1);

      if (fabs (f_i) < GSL_DBL_EPSILON * (fabs (d_i) + fabs (d_ip1)))
        {
          gsl_vector_set (f, i, 0.0);
        }
      d_i = d_ip1;
    }

}

void qrstep (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V)
{ 
#if !USE_BLAS
  const size_t M = U->size1;
  const size_t N = V->size1;
#endif
  const size_t n = d->size;
  double y, z;
  double ak, bk, zk, ap, bp, aq, bq;
  size_t i, k;

  if (n == 1)
    return;  /* shouldn't happen */

  /* Compute 2x2 svd directly */

  if (n == 2)
    {
      svd2 (d, f, U, V);
      return;
    }

  /* Chase out any zeroes on the diagonal */

  for (i = 0; i < n - 1; i++)
    {
      double d_i = gsl_vector_get (d, i);
      
      if (d_i == 0.0)
        {
          chase_out_intermediate_zero (d, f, U, i);
          return;
        }
    }

  /* Chase out any zero at the end of the diagonal */

  {
    double d_nm1 = gsl_vector_get (d, n - 1);

    if (d_nm1 == 0.0) 
      {
        chase_out_trailing_zero (d, f, V);
        return;
      }
  }


  /* Apply QR reduction steps to the diagonal and offdiagonal */

  {
    double d0 = gsl_vector_get (d, 0);
    double f0 = gsl_vector_get (f, 0);
    
    double d1 = gsl_vector_get (d, 1);
    double f1 = gsl_vector_get (f, 1);
    
    {
      double mu = trailing_eigenvalue (d, f);
    
      y = d0 * d0 - mu;
      z = d0 * f0;
    }
    
    /* Set up the recurrence for Givens rotations on a bidiagonal matrix */
    
    ak = 0;
    bk = 0;
    
    ap = d0;
    bp = f0;
    
    aq = d1;
    bq = f1;
  }

  for (k = 0; k < n - 1; k++)
    {
      double c, s;
      create_givens (y, z, &c, &s);

      /* Compute V <= V G */

#ifdef USE_BLAS
      {
        gsl_vector_view Vk = gsl_matrix_column(V,k);
        gsl_vector_view Vkp1 = gsl_matrix_column(V,k+1);
        gsl_blas_drot(&Vk.vector, &Vkp1.vector, c, -s);
      }
#else
      for (i = 0; i < N; i++)
        {
          double Vip = gsl_matrix_get (V, i, k);
          double Viq = gsl_matrix_get (V, i, k + 1);
          gsl_matrix_set (V, i, k, c * Vip - s * Viq);
          gsl_matrix_set (V, i, k + 1, s * Vip + c * Viq);
        }
#endif

      /* compute B <= B G */

      {
        double bk1 = c * bk - s * z;

        double ap1 = c * ap - s * bp;
        double bp1 = s * ap + c * bp;
        double zp1 = -s * aq;

        double aq1 = c * aq;

        if (k > 0)
          {
            gsl_vector_set (f, k - 1, bk1);
          }

        ak = ap1;
        bk = bp1;
        zk = zp1;

        ap = aq1;

        if (k < n - 2)
          {
            bp = gsl_vector_get (f, k + 1);
          }
        else
          {
            bp = 0.0;
          }

        y = ak;
        z = zk;
      }

      create_givens (y, z, &c, &s);

      /* Compute U <= U G */

#ifdef USE_BLAS
      {
        gsl_vector_view Uk = gsl_matrix_column(U,k);
        gsl_vector_view Ukp1 = gsl_matrix_column(U,k+1);
        gsl_blas_drot(&Uk.vector, &Ukp1.vector, c, -s);
      }
#else
      for (i = 0; i < M; i++)
        {
          double Uip = gsl_matrix_get (U, i, k);
          double Uiq = gsl_matrix_get (U, i, k + 1);
          gsl_matrix_set (U, i, k, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
        }
#endif

      /* compute B <= G^T B */

      {
        double ak1 = c * ak - s * zk;
        double bk1 = c * bk - s * ap;
        double zk1 = -s * bp;

        double ap1 = s * bk + c * ap;
        double bp1 = c * bp;

        gsl_vector_set (d, k, ak1);

        ak = ak1;
        bk = bk1;
        zk = zk1;

        ap = ap1;
        bp = bp1;

        if (k < n - 2)
          {
            aq = gsl_vector_get (d, k + 2);
          }
        else
          {
            aq = 0.0;
          }

        y = bk;
        z = zk;
      }
    }

  gsl_vector_set (f, n - 2, bk);
  gsl_vector_set (d, n - 1, ap);
}

double trailing_eigenvalue (const gsl_vector * d, const gsl_vector * f)
{

  const size_t n = d->size;

  double da = gsl_vector_get (d, n - 2);
  double db = gsl_vector_get (d, n - 1);
  double fa = (n > 2) ? gsl_vector_get (f, n - 3) : 0.0;
  double fb = gsl_vector_get (f, n - 2);

  double ta = da * da + fa * fa;
  double tb = db * db + fb * fb;
  double tab = da * fb;

  double dt = (ta - tb) / 2.0;

  double mu;

  if (dt >= 0)
    {
      mu = tb - (tab * tab) / (dt + hypot (dt, tab));
    }
  else 
    {
      mu = tb + (tab * tab) / ((-dt) + hypot (dt, tab));
    }

  return mu;
}

void create_schur (double d0, double f0, double d1, double * c, double * s)
{

  double apq = 2.0 * d0 * f0;
  
  if (apq != 0.0)
    {
      double t;
      double tau = (f0*f0 + (d1 + d0)*(d1 - d0)) / apq;
      
      if (tau >= 0.0)
        {
          t = 1.0/(tau + hypot(1.0, tau));
        }
      else
        {
          t = -1.0/(-tau + hypot(1.0, tau));
        }

      *c = 1.0 / hypot(1.0, t);
      *s = t * (*c);
    }
  else
    {
      *c = 1.0;
      *s = 0.0;
    }
}

void svd2 (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V)
{

  size_t i;
  double c, s, a11, a12, a21, a22;

  const size_t M = U->size1;
  const size_t N = V->size1;

  double d0 = gsl_vector_get (d, 0);
  double f0 = gsl_vector_get (f, 0);
  
  double d1 = gsl_vector_get (d, 1);

  if (d0 == 0.0)
    {
      /* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */

      create_givens (f0, d1, &c, &s);

      /* compute B <= G^T B X,  where X = [0,1;1,0] */

      gsl_vector_set (d, 0, c * f0 - s * d1);
      gsl_vector_set (f, 0, s * f0 + c * d1);
      gsl_vector_set (d, 1, 0.0);
      
      /* Compute U <= U G */

      for (i = 0; i < M; i++)
        {
          double Uip = gsl_matrix_get (U, i, 0);
          double Uiq = gsl_matrix_get (U, i, 1);
          gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
        }

      /* Compute V <= V X */

      gsl_matrix_swap_columns (V, 0, 1);

      return;
    }
  else if (d1 == 0.0)
    {
      /* Eliminate off-diagonal element in [d0,f0;0,0] */

      create_givens (d0, f0, &c, &s);

      /* compute B <= B G */

      gsl_vector_set (d, 0, d0 * c - f0 * s);
      gsl_vector_set (f, 0, 0.0);

      /* Compute V <= V G */

      for (i = 0; i < N; i++)
        {
          double Vip = gsl_matrix_get (V, i, 0);
          double Viq = gsl_matrix_get (V, i, 1);
          gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
          gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
        }

      return;
    }
  else
    {
      /* Make columns orthogonal, A = [d0, f0; 0, d1] * G */
      
      create_schur (d0, f0, d1, &c, &s);
      
      /* compute B <= B G */
      
      a11 = c * d0 - s * f0;
      a21 = - s * d1;
      
      a12 = s * d0 + c * f0;
      a22 = c * d1;
      
      /* Compute V <= V G */
      
      for (i = 0; i < N; i++)
        {
          double Vip = gsl_matrix_get (V, i, 0);
          double Viq = gsl_matrix_get (V, i, 1);
          gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
          gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
        }
      
      /* Eliminate off-diagonal elements, bring column with largest
         norm to first column */
      
      if (hypot(a11, a21) < hypot(a12,a22))
        {
          double t1, t2;

          /* B <= B X */

          t1 = a11; a11 = a12; a12 = t1;
          t2 = a21; a21 = a22; a22 = t2;

          /* V <= V X */

          gsl_matrix_swap_columns(V, 0, 1);
        } 

      create_givens (a11, a21, &c, &s);
      
      /* compute B <= G^T B */
      
      gsl_vector_set (d, 0, c * a11 - s * a21);
      gsl_vector_set (f, 0, c * a12 - s * a22);
      gsl_vector_set (d, 1, s * a12 + c * a22);
      
      /* Compute U <= U G */
      
      for (i = 0; i < M; i++)
        {
          double Uip = gsl_matrix_get (U, i, 0);
          double Uiq = gsl_matrix_get (U, i, 1);
          gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
          gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
        }

      return;
    }
}


void chase_out_intermediate_zero (gsl_vector * d, gsl_vector * f, gsl_matrix * U, size_t k0)
{
  
#if !USE_BLAS
  const size_t M = U->size1;
#endif
  const size_t n = d->size;
  double c, s;
  double x, y;
  size_t k;

  x = gsl_vector_get (f, k0);
  y = gsl_vector_get (d, k0+1);

  for (k = k0; k < n - 1; k++)
    {
      create_givens (y, -x, &c, &s);
      
      /* Compute U <= U G */

#ifdef USE_BLAS
      {
        gsl_vector_view Uk0 = gsl_matrix_column(U,k0);
        gsl_vector_view Ukp1 = gsl_matrix_column(U,k+1);
        gsl_blas_drot(&Uk0.vector, &Ukp1.vector, c, -s);
      }
#else
      {
        size_t i;

        for (i = 0; i < M; i++)
          {
            double Uip = gsl_matrix_get (U, i, k0);
            double Uiq = gsl_matrix_get (U, i, k + 1);
            gsl_matrix_set (U, i, k0, c * Uip - s * Uiq);
            gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
          }
      }
#endif
      
      /* compute B <= G^T B */
      
      gsl_vector_set (d, k + 1, s * x + c * y);

      if (k == k0)
        gsl_vector_set (f, k, c * x - s * y );

      if (k < n - 2) 
        {
          double z = gsl_vector_get (f, k + 1);
          gsl_vector_set (f, k + 1, c * z); 

          x = -s * z ;
          y = gsl_vector_get (d, k + 2); 
        }
    }
}

void chase_out_trailing_zero (gsl_vector * d, gsl_vector * f, gsl_matrix * V)
{
#if !USE_BLAS
  const size_t N = V->size1;
#endif
  const size_t n = d->size;
  double c, s;
  double x, y;
  size_t k;

  x = gsl_vector_get (d, n - 2);
  y = gsl_vector_get (f, n - 2);

  for (k = n - 1; k > 0 && k--;)
    {
      create_givens (x, y, &c, &s);

      /* Compute V <= V G where G = [c, s ; -s, c] */

#ifdef USE_BLAS
      {
        gsl_vector_view Vp = gsl_matrix_column(V,k);
        gsl_vector_view Vq = gsl_matrix_column(V,n-1);
        gsl_blas_drot(&Vp.vector, &Vq.vector, c, -s);
      }
#else
      {
        size_t i;
   
        for (i = 0; i < N; i++)
          {
            double Vip = gsl_matrix_get (V, i, k);
            double Viq = gsl_matrix_get (V, i, n - 1);
            gsl_matrix_set (V, i, k, c * Vip - s * Viq);
            gsl_matrix_set (V, i, n - 1, s * Vip + c * Viq);
          }
      }
#endif

      /* compute B <= B G */
      
      gsl_vector_set (d, k, c * x - s * y);

      if (k == n - 2)
        gsl_vector_set (f, k, s * x + c * y );

      if (k > 0) 
        {
          double z = gsl_vector_get (f, k - 1);
          gsl_vector_set (f, k - 1, c * z); 

          x = gsl_vector_get (d, k - 1); 
          y = s * z ;
        }
    }
}

int gsl_linalg_bidiag_unpack (const gsl_matrix * A, const gsl_vector * tau_U, 
			      gsl_matrix * U, const gsl_vector * tau_V,
			      gsl_matrix * V, gsl_vector * diag, 
			      gsl_vector * superdiag) {
  const size_t M = A->size1;
  const size_t N = A->size2;

  const size_t K = GSL_MIN(M, N);

  if (M < N)
    {
      bpm_error("matrix A must have M >= N in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (tau_U->size != K)
    {
      bpm_error("size of tau must be MIN(M,N) in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (tau_V->size + 1 != K)
    {
      bpm_error("size of tau must be MIN(M,N) - 1 in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (U->size1 != M || U->size2 != N)
    {
      bpm_error("size of U must be M x N in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (V->size1 != N || V->size2 != N)
    {
      bpm_error("size of V must be N x N in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (diag->size != K)
    {
      bpm_error("size of diagonal must match size of A in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (superdiag->size + 1 != K)
    {
      bpm_error("size of subdiagonal must be (diagonal size - 1) in gsl_linalg_bidiag_unpack(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else
    {
      size_t i, j;

      /* Copy diagonal into diag */

      for (i = 0; i < N; i++)
        {
          double Aii = gsl_matrix_get (A, i, i);
          gsl_vector_set (diag, i, Aii);
        }

      /* Copy superdiagonal into superdiag */

      for (i = 0; i < N - 1; i++)
        {
          double Aij = gsl_matrix_get (A, i, i+1);
          gsl_vector_set (superdiag, i, Aij);
        }

      /* Initialize V to the identity */

      gsl_matrix_set_identity (V);

      for (i = N - 1; i > 0 && i--;)
        {
          /* Householder row transformation to accumulate V */
          gsl_vector_const_view r = gsl_matrix_const_row (A, i);
          gsl_vector_const_view h = 
            gsl_vector_const_subvector (&r.vector, i + 1, N - (i+1));
          
          double ti = gsl_vector_get (tau_V, i);
          
          gsl_matrix_view m = 
            gsl_matrix_submatrix (V, i + 1, i + 1, N-(i+1), N-(i+1));
          
          gsl_linalg_householder_hm (ti, &h.vector, &m.matrix);
        }

      /* Initialize U to the identity */

      gsl_matrix_set_identity (U);

      for (j = N; j > 0 && j--;)
        {
          /* Householder column transformation to accumulate U */
          gsl_vector_const_view c = gsl_matrix_const_column (A, j);
          gsl_vector_const_view h = gsl_vector_const_subvector (&c.vector, j, M - j);
          double tj = gsl_vector_get (tau_U, j);
          
          gsl_matrix_view m = 
            gsl_matrix_submatrix (U, j, j, M-j, N-j);
          
          gsl_linalg_householder_hm (tj, &h.vector, &m.matrix);
        }

      return BPM_SUCCESS;
    }
}

int gsl_linalg_bidiag_unpack2 (gsl_matrix * A, gsl_vector * tau_U, gsl_vector * tau_V, gsl_matrix * V)
{

  const size_t M = A->size1;
  const size_t N = A->size2;

  const size_t K = GSL_MIN(M, N);

  if (M < N)
    {
      bpm_error("matrix A must have M >= N in gsl_linalg_bidiag_unpack2(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (tau_U->size != K)
    {
      bpm_error("size of tau must be MIN(M,N) in gsl_linalg_bidiag_unpack2(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (tau_V->size + 1 != K)
    {
      bpm_error("size of tau must be MIN(M,N) - 1 in gsl_linalg_bidiag_unpack2(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else if (V->size1 != N || V->size2 != N)
    {
      bpm_error("size of V must be N x N in gsl_linalg_bidiag_unpack2(...)",
		__FILE__, __LINE__);
      return BPM_FAILURE;
    }
  else
    {
      size_t i, j;

      /* Initialize V to the identity */

      gsl_matrix_set_identity (V);

      for (i = N - 1; i > 0 && i--;)
        {
          /* Householder row transformation to accumulate V */
          gsl_vector_const_view r = gsl_matrix_const_row (A, i);
          gsl_vector_const_view h = 
            gsl_vector_const_subvector (&r.vector, i + 1, N - (i+1));
          
          double ti = gsl_vector_get (tau_V, i);
          
          gsl_matrix_view m = 
            gsl_matrix_submatrix (V, i + 1, i + 1, N-(i+1), N-(i+1));
          
          gsl_linalg_householder_hm (ti, &h.vector, &m.matrix);
        }

      /* Copy superdiagonal into tau_v */

      for (i = 0; i < N - 1; i++)
        {
          double Aij = gsl_matrix_get (A, i, i+1);
          gsl_vector_set (tau_V, i, Aij);
        }

      /* Allow U to be unpacked into the same memory as A, copy
         diagonal into tau_U */

      for (j = N; j > 0 && j--;)
        {
          /* Householder column transformation to accumulate U */
          double tj = gsl_vector_get (tau_U, j);
          double Ajj = gsl_matrix_get (A, j, j);
          gsl_matrix_view m = gsl_matrix_submatrix (A, j, j, M-j, N-j);

          gsl_vector_set (tau_U, j, Ajj);
          gsl_linalg_householder_hm1 (tj, &m.matrix);
        }

      return BPM_SUCCESS;
    }
}

int gsl_linalg_SV_decomp (gsl_matrix * A, gsl_matrix * V, gsl_vector * S, 
			  gsl_vector * work) {

  size_t a, b, i, j;

  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t K = GSL_MIN (M, N);

  if (M < N)
    {
      bpm_error( "svd of MxN matrix, M<N, is not implemented in gsl_linalg_SV_solve(...)",
		 __FILE__, __LINE__);
    }
  else if (V->size1 != N)
    {
      bpm_error( "square matrix V must match second dimension of matrix A in gsl_linalg_SV_solve(...)",
		 __FILE__, __LINE__);
    }
  else if (V->size1 != V->size2)
    {
      bpm_error( "matrix V must be square in gsl_linalg_SV_solve(...)",
		 __FILE__, __LINE__);
    }
  else if (S->size != N)
    {
      bpm_error( "length of vector S must match second dimension of matrix A in gsl_linalg_SV_solve(...)",
		 __FILE__, __LINE__);
    }
  else if (work->size != N)
    {
      bpm_error( "length of workspace must match second dimension of matrix A in gsl_linalg_SV_solve(...)",
		 __FILE__, __LINE__);
    }

  /* Handle the case of N = 1 (SVD of a column vector) */

  if (N == 1)
    {
      gsl_vector_view column = gsl_matrix_column (A, 0);
      double norm = gsl_blas_dnrm2 (&column.vector);

      gsl_vector_set (S, 0, norm); 
      gsl_matrix_set (V, 0, 0, 1.0);
      
      if (norm != 0.0)
        {
          gsl_blas_dscal (1.0/norm, &column.vector);
        }

      return BPM_SUCCESS;
    }
  
  {
    gsl_vector_view f = gsl_vector_subvector (work, 0, K - 1);
    
    /* bidiagonalize matrix A, unpack A into U S V */

    gsl_linalg_bidiag_decomp (A, S, &f.vector);
    gsl_linalg_bidiag_unpack2 (A, S, &f.vector, V);

    /* apply reduction steps to B=(S,Sd) */
    chop_small_elements (S, &f.vector);

    /* Progressively reduce the matrix until it is diagonal */
    
    b = N - 1;
    
    while (b > 0)
      {
        double fbm1 = gsl_vector_get (&f.vector, b - 1);

        if (fbm1 == 0.0 || gsl_isnan (fbm1))
          {
            b--;
            continue;
          }

        /* Find the largest unreduced block (a,b) starting from b
           and working backwards */
        
        a = b - 1;
        
        while (a > 0)
          {
            double fam1 = gsl_vector_get (&f.vector, a - 1);

            if (fam1 == 0.0 || gsl_isnan (fam1))
              {
                break;
              }
            
            a--;

          }

        {
          const size_t n_block = b - a + 1;

          gsl_vector_view S_block = gsl_vector_subvector (S, a, n_block);
          gsl_vector_view f_block = gsl_vector_subvector (&f.vector, a, n_block - 1);
          
          gsl_matrix_view U_block =
            gsl_matrix_submatrix (A, 0, a, A->size1, n_block);
          gsl_matrix_view V_block =
            gsl_matrix_submatrix (V, 0, a, V->size1, n_block);
          
          qrstep (&S_block.vector, &f_block.vector, &U_block.matrix, &V_block.matrix);
          /* remove any small off-diagonal elements */
          
          chop_small_elements (&S_block.vector, &f_block.vector);
        }
      }
  }
  /* Make singular values positive by reflections if necessary */
  
  for (j = 0; j < K; j++)
    {
      double Sj = gsl_vector_get (S, j);
      
      if (Sj < 0.0)
        {
          for (i = 0; i < N; i++)
            {
              double Vij = gsl_matrix_get (V, i, j);

              gsl_matrix_set (V, i, j, -Vij);
            }
          
          gsl_vector_set (S, j, -Sj);
        }
    }
  
  /* Sort singular values into decreasing order */
  
  for (i = 0; i < K; i++)
    {
      double S_max = gsl_vector_get (S, i);
      size_t i_max = i;
      
      for (j = i + 1; j < K; j++)
        {
          double Sj = gsl_vector_get (S, j);
          
          if (Sj > S_max)
            {
              S_max = Sj;
              i_max = j;
            }
        }
      
      if (i_max != i)
        {
          /* swap eigenvalues */
          gsl_vector_swap_elements (S, i, i_max);
          
          /* swap eigenvectors */
          gsl_matrix_swap_columns (A, i, i_max);
          gsl_matrix_swap_columns (V, i, i_max);
        }
    }
  
  return BPM_SUCCESS;
}

