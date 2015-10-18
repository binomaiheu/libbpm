/**
  @file
  @ingroup message
*/
#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_interface.h>

void bpm_warning( char* msg, char* f, int l ) {

  fprintf( stderr, "+++ libbpm warning (%s <%d>) [evt:%d] :: %s \n",f, l, libbpm_evtnum, msg );
  fflush( stderr );

  return;
}
