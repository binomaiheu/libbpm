/**
   @file
   @ingroup analysis
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_analysis.h>

int (*ana_cutfn)( bpmproc_t *proc ) = ana_def_cutfn;

int ana_def_cutfn( bpmproc_t *proc ) {

  return BPM_GOOD_EVENT;
}
