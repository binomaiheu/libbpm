/** 
    @file     
    @ingroup message
    @brief   libbpm error/warning messages

    This header defines the routines which take care of printing 
    error and warning messages
*/

/** @addtogroup message */    
/** @{ */

#ifndef BPMERR_H__
#define BPMERR_H__

/* -----------------------------------------------------------------------------
// includes
// -------------------------------------------------------------------------- */
#include <bpm/bpm_defs.h>


/* -----------------------------------------------------------------------------
// macro definitions
// -------------------------------------------------------------------------- */


/* -----------------------------------------------------------------------------
// typedefs, enums and other declarations
// -------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif


/* -----------------------------------------------------------------------------
// function prototypes and declarations
// -------------------------------------------------------------------------- */

  /**
    Prints an error message in a standard format

    @param msg the error messages, without end of line character
    @param f the file position (__FILE__)
    @param l the line in the file (__LINE__)

    @return void
  */
  EXTERN void bpm_error( char* msg, char* f, int l );


  /**
    Prints an warning message in a standard format
    
    @param msg the error messages, without end of line character
    @param f the file position (__FILE__)
    @param l the line in the file (__LINE__)
    
    @return void
  */
  EXTERN void bpm_warning( char* msg, char* f, int l );


#ifdef __cplusplus
}
#endif

#endif /* #ifndef BPMERR_H__ */
/*@}*/
/* ================================ end of file ============================= */
