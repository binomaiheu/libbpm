/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Swap two matrix columns

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

    @param m The matrix
    @param i index of column one
    @param j index of column two
    
    @return BPM_SUCCESS if everything was OK, BPM_FAILURE if not
*/

int gsl_matrix_swap_columns(gsl_matrix * m, const size_t i, const size_t j) {
  
  const size_t size1 = m->size1;
  const size_t size2 = m->size2;

  if (i >= size2)
    {
      bpm_error("first column index is out of range in gsl_matrix_swap_columns(...)",
		__FILE__, __LINE__ );
      return BPM_FAILURE;
    }

  if (j >= size2)
    {
      bpm_error("second column index is out of range in gsl_matrix_swap_columns(...)",
		__FILE__, __LINE__ );
      return BPM_FAILURE;
    }

  if (i != j)
    {
      double *col1 = m->data + i;
      double *col2 = m->data + j;
      
      size_t p;
      
      for (p = 0; p < size1; p++)
        {
          size_t k;
          size_t n = p * m->tda;
 
          for (k = 0; k < 1; k++)
            {
              double tmp = col1[n+k] ;
              col1[n+k] = col2[n+k] ;
              col2[n+k] = tmp ;
            }
        }
    }

  return BPM_SUCCESS;
}


// -------------------------------------------------------------------

/**
   Retrieve a column of a matrix

   @param m The matrix
   @param j index of the column
   
   @return BPM_SUCCESS if everything was OK, BPM_FAILURE if not
*/

_gsl_vector_view gsl_matrix_column (gsl_matrix * m, const size_t j) {
 
  _gsl_vector_view view = NULL_VECTOR_VIEW;
  
  if (j >= m->size2)
    {
      bpm_error("column index is out of range in gsl_matrix_column", 
		__FILE__, __LINE__ );
      return view;
    }

  {
    gsl_vector v = NULL_VECTOR;
    
    v.data = m->data + j;
    v.size = m->size1;
    v.stride = m->tda;
    v.block = m->block;
    v.owner = 0;

    view.vector = v;
    return view;
  }
}

/**
   Get the matrix value associated with the given row and column
   
   @param m The matrix
   @param i The row number
   @param j The column number
   
   @return The value of the matrix element
*/
double gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j) {

  return *(double*) (m->data + (i * m->tda + j));
}




/**
   Set the matrix value associated with the given row and column
   
   @param m The matrix
   @param i The row number
   @param j The column number
   @param x the value to set
   
*/
void gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x) { 

  *(double*) (m->data + (i * m->tda + j)) = x;

}



/**
   Retrieve a submatrix of the given matrix
*/
_gsl_matrix_view gsl_matrix_submatrix (gsl_matrix * m, const size_t i, const size_t j, 
				       const size_t n1, const size_t n2) {

  _gsl_matrix_view view = NULL_MATRIX_VIEW; 

  if (i >= m->size1)
    {
      bpm_error("row index is out of range in gsl_matrix_submatrix(...)",
		__FILE__, __LINE__ );
      return view;
    }
  else if (j >= m->size2)
    {
      bpm_error("column index is out of range in gsl_matrix_submatrix(...)",
		__FILE__, __LINE__ );
      return view;
    }
  else if (n1 == 0)
    {
      bpm_error("first dimension must be non-zero in gsl_matrix_submatrix(...)",
		__FILE__, __LINE__ );
      return view;
    }
  else if (n2 == 0)
    {
      bpm_error("second dimension must be non-zero in gsl_matrix_submatrix(...)",
		__FILE__, __LINE__ );
      return view;
    }
  else if (i + n1 > m->size1)
    {
      bpm_error("first dimension overflows matrix in gsl_matrix_submatrix(...)",
		__FILE__, __LINE__ );
      return view;
    }
  else if (j + n2 > m->size2)
    {
      bpm_error("second dimension overflows matrix in gsl_matrix_submatrix(...)",
		__FILE__, __LINE__ );
      return view;
    }

  {
     gsl_matrix s = NULL_MATRIX;

     s.data = m->data + (i * m->tda + j);
     s.size1 = n1;
     s.size2 = n2;
     s.tda = m->tda;
     s.block = m->block;
     s.owner = 0;
     
     view.matrix = s;     
     return view;
  }
}

gsl_matrix * gsl_matrix_alloc (const size_t n1, const size_t n2)
{

  gsl_block * block;
  gsl_matrix * m;

  if (n1 == 0)
    {
      bpm_error("matrix dimension n1 must be positive integer in gsl_matrix_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }
  else if (n2 == 0)
    {
      bpm_error("matrix dimension n2 must be positive integer in gsl_matrix_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }

  m = (gsl_matrix *) malloc (sizeof(gsl_matrix));

  if (m == 0)
    {
      bpm_error("failed to allocate space for matrix struct in gsl_matrix_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }

  /* FIXME: n1*n2 could overflow for large dimensions */

  block = gsl_block_alloc (n1 * n2) ;

  if (block == 0)
    {
      bpm_error("failed to allocate space for block in gsl_matrix_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }

  m->data = block->data;
  m->size1 = n1;
  m->size2 = n2;
  m->tda = n2; 
  m->block = block;
  m->owner = 1;

  return m;
}

gsl_matrix * gsl_matrix_calloc (const size_t n1, const size_t n2)
{
  size_t i;

  gsl_matrix * m = gsl_matrix_alloc( n1, n2);

  if (m == 0)
    return 0;

  /* initialize matrix to zero */

  for (i = 0; i < n1 * n2; i++)
    {
      m->data[i] = 0;
    }

  return m;
}

_gsl_vector_const_view gsl_matrix_const_row (const gsl_matrix * m, const size_t i)
{

  _gsl_vector_const_view view = NULL_VECTOR_VIEW;
  
  if (i >= m->size1)
    {
      bpm_error("row index is out of range in gsl_matrix_const_row(...)",
		__FILE__, __LINE__);
      return view;
    }
  
  {
    gsl_vector v = NULL_VECTOR;
    
    v.data = m->data + i * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;
    
    view.vector = v;
    return view;
  }
}

_gsl_vector_view gsl_matrix_row (gsl_matrix * m, const size_t i)
{

  _gsl_vector_view view = NULL_VECTOR_VIEW;
  
  if (i >= m->size1)
    {
      bpm_error("row index is out of range in gsl_matrix_row(...)",
		__FILE__, __LINE__);
      return view;
    }
  
  {
    gsl_vector v = NULL_VECTOR;
    
    v.data = m->data + i * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;
    
    view.vector = v;
    return view;
  }
}

_gsl_vector_const_view gsl_matrix_const_column (const gsl_matrix * m, const size_t i)
{

  _gsl_vector_const_view view = NULL_VECTOR_VIEW;
  
  if (i >= m->size1)
    {
      bpm_error("row index is out of range", __FILE__, __LINE__ );
      return view;
    }
  
  {
    gsl_vector v = NULL_VECTOR;
    
    v.data = m->data + i * m->tda;
    v.size = m->size2;
    v.stride = 1;
    v.block = m->block;
    v.owner = 0;
    
    view.vector = v;
    return view;
  }
}

void gsl_matrix_set_identity (gsl_matrix * m)
{

  size_t i, j;
  double * const data = m->data;
  const size_t p = m->size1 ;
  const size_t q = m->size2 ;
  const size_t tda = m->tda ;

  const double zero = 0.;
  const double one = 1.;

  for (i = 0; i < p; i++)
    {
      for (j = 0; j < q; j++)
        {
          *(double *) (data + (i * tda + j)) = ((i == j) ? 1. : 0.);
        }
    }
}

