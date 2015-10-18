/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

gsl_block *gsl_block_alloc(const size_t n)
{
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
  
  gsl_block * b;

  if (n == 0)
    {
      bpm_error("block length n must be positive integer in gsl_block_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }

  b = (gsl_block *) malloc (sizeof (gsl_block));

  if (b == 0)
    {
      bpm_error("failed to allocate space for block struct in gsl_block_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }

  b->data = (double *) malloc (n * sizeof (double));

  if (b->data == 0)
    {
      free (b);         /* exception in constructor, avoid memory leak */
      
      bpm_error("failed to allocate space for block data in gsl_block_alloc(...)",
		__FILE__, __LINE__);
      return NULL;
    }

  b->size = n;

  return b;
}

void gsl_block_free (gsl_block * b)
{

  free (b->data);
  free (b);
}
