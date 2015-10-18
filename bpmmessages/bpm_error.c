/**
  @file
  @ingroup message
*/
#include <stdio.h>
#include <bpm/bpm_messages.h>
#include <bpm/bpm_interface.h>

void bpm_error( char* msg, char* f, int l ) {

  fprintf( stderr, "*** libbpm error (%s <%d>) [evt:%d] :: %s \n",f, l, libbpm_evtnum, msg );
  fflush( stderr );

  return;
}
