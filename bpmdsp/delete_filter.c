/**
   @file
   @ingroup dsp
*/
#include "bpm/bpm_dsp.h"

void delete_filter( filter_t *f ) {
  if(!f)return;
  if(f->cplane) free(f->cplane);
  if(f->wfbuffer) free(f->wfbuffer );
  free(f);
  return;
}
