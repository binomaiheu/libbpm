/**
   @file
   @ingroup analysis
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_analysis.h>

int ana_set_cutfn( int (*cutfn)( bpmproc_t *proc ) ) {
  
  if (! cutfn ){
    
    bpm_error("Bad function pointer in ana_set_ctufn(...)",
	      __FILE__, __LINE__);
  }

  ana_cutfn = cutfn;

  return BPM_GOOD_EVENT;
}
