/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

_gsl_vector_view gsl_vector_subvector (gsl_vector *v, size_t offset, size_t n) {

  /**
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

  _gsl_vector_view view = NULL_VECTOR_VIEW;

  if (n == 0)
    {
      bpm_error("vector length n must be positive integer in gsl_vector_subvector(...)",
		__FILE__, __LINE__);
      return view;
    }

  if (offset + (n - 1) >= v->size)
    {
      bpm_error("view would extend past end of vector in gsl_vector_subvector(...)", 
		__FILE__, __LINE__);
      return view;
    }

  {
    gsl_vector s = NULL_VECTOR;

    s.data = v->data +  v->stride * offset ;
    s.size = n;
    s.stride = v->stride;
    s.block = v->block;
    s.owner = 0;

    view.vector = s;
    return view;
  }
}

double gsl_vector_get (const gsl_vector * v, const size_t i)
{
  /**
    The following line is a generalization of return v->data[i]
  */

  return *(double *) (v->data + i * v->stride);
}

void gsl_vector_set (gsl_vector * v, const size_t i, double x)
{
  /**
     The following line is a generalization of v->data[i] = x 
  */

  *(double *) (v->data + i * v->stride) = x;
}

int gsl_vector_swap_elements (gsl_vector * v, const size_t i, const size_t j) {
  

  double * data = v->data ;
  const size_t size = v->size ;
  const size_t stride = v->stride ;

  if (i >= size)
    {
      bpm_error("first index is out of range in gsl_vector_swap_elements(...)",
		__FILE__, __LINE__ );
      return BPM_FAILURE;
    }

  if (j >= size)
    {
            bpm_error("second index is out of range in gsl_vector_swap_elements(...)",
		__FILE__, __LINE__ );
      return BPM_FAILURE;
    }

  if (i != j)
    {
      const size_t s = stride ;
      size_t k ;

      for (k = 0; k < 1; k++)
        {
          double tmp = data[j*s + k];
          data[j*s+k] = data[i*s + k];
          data[i*s+k] = tmp;
        }
    }
  
  return BPM_SUCCESS;
}

gsl_vector *gsl_vector_alloc (const size_t n)
{
 
  gsl_block * block;
  gsl_vector * v;

  if (n == 0)
    {
      bpm_error("vector length n must be positive integer in gsl_vector_alloc(...)",
		__FILE__, __LINE__ );
      return NULL;
    }

  v = (gsl_vector*) malloc (sizeof (gsl_vector));

  if (v == 0)
    {
      bpm_error("failed to allocate space for vector struct in gsl_vector_alloc(...)",
		__FILE__, __LINE__ );
      return NULL;
    }

  block = gsl_block_alloc(n);

  if (block == 0)
    {
      free (v) ;
      bpm_error("failed to allocate space for block in gsl_vector_alloc(...)",
		__FILE__, __LINE__ );
      return NULL;
    }
      
  v->data = block->data ;
  v->size = n;
  v->stride = 1;
  v->block = block;
  v->owner = 1;

  return v;
}

gsl_vector *gsl_vector_calloc (const size_t n)
{
  
  size_t i;

  gsl_vector* v = gsl_vector_alloc (n);

  if (v == 0)
    return 0;

  /* initialize vector to zero */

  for (i = 0; i < n; i++)
    {
      v->data[i] = 0;
    }

  return v;
}


_gsl_vector_const_view gsl_vector_const_subvector (const gsl_vector *v, size_t offset, size_t n)
{
  
  _gsl_vector_const_view view = NULL_VECTOR_VIEW;

  if (n == 0)
    {
      bpm_error("vector length n must be positive integer in gsl_vector_const_subvector(...)",
		__FILE__, __LINE__);
      return view;
    }

  if (offset + (n - 1) >= v->size)
    {
      bpm_error("view would extend past end of vector in gsl_vector_const_subvector(...)", 
		__FILE__, __LINE__);
      return view;
    }

  {
    gsl_vector s = NULL_VECTOR;

    s.data = v->data +  v->stride * offset ;
    s.size = n;
    s.stride = v->stride;
    s.block = v->block;
    s.owner = 0;

    view.vector = s;
    return view;
  }
}

void gsl_vector_free (gsl_vector * v)
{

  if (v->owner)
    {
      gsl_block_free(v->block) ;
    }
  free (v);
}
