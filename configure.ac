dnl Process this file with autoconf to produce a configure script.

dnl Library versioning : BPM_LIB_VERSIONING
dnl distribution  so-version
dnl libbpm-0.1    0.1  (CURRENT:REVISION)
dnl libbpm-0.2    0.2
dnl libbpm-0.3    0.3
AC_INIT([libbpm],[0.3])
AC_CONFIG_SRCDIR(bpm_defs.h)
BPM_LIB_VERSION="0.3"
AC_SUBST(BPM_LIB_VERSION)

AM_INIT_AUTOMAKE([gnu no-dependencies])
dnl things required by automake
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_LANG(C)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S

dnl Disable unnecessary libtool tests for c++,fortran,java
define([AC_LIBTOOL_LANG_CXX_CONFIG], [:])dnl
define([AC_LIBTOOL_LANG_F77_CONFIG], [:])dnl
define([AC_LIBTOOL_LANG_GCJ_CONFIG], [:])dnl

dnl libtool:
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

dnl Check compiler features
AC_TYPE_SIZE_T
AC_C_CONST

BPM_CFLAGS="-I$includedir"
BPM_LIBS="-L$libdir -lbpm -lbpmprocess -lbpmanalysis -lbpmcalibration"
BPM_LIBS="$BPM_LIBS -lbpminterface -lbpmmessages -lbpmnr -lbpmorbit -lbpmrf -lbpmsimulation"
BPM_LIBS="$BPM_LIBS -lbpmdsp -lbpmwf"

AC_SUBST(BPM_CFLAGS)
AC_SUBST(BPM_LIBS)

dnl Checks for header files.

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.

dnl The doxygen stuff
dnl ---------------------------------------------------------------------
DX_HTML_FEATURE(ON)
DX_PS_FEATURE(ON)
DX_PDF_FEATURE(ON)
DX_INIT_DOXYGEN(libbpm, doxygen.cfg, documentation)

AM_CONDITIONAL([DX_COND_doc],  false)
AM_CONDITIONAL([DX_COND_dot],  false)
AM_CONDITIONAL([DX_COND_html], true)
AM_CONDITIONAL([DX_COND_ps],   true)
AM_CONDITIONAL([DX_COND_pdf],  true)
AM_CONDITIONAL([DX_COND_chm],  false)
AM_CONDITIONAL([DX_COND_chi],  false)
AM_CONDITIONAL([DX_COND_man],  false)
AM_CONDITIONAL([DX_COND_rtf],  false)
AM_CONDITIONAL([DX_COND_xml],  false)

dnl
AC_CONFIG_FILES([libbpm-config libbpm.pc bpm_version.h Makefile \
                bpm/Makefile \
                bpmanalysis/Makefile \
		bpmcalibration/Makefile \
                bpminterface/Makefile \
                bpmmessages/Makefile \
                bpmnr/Makefile \
                bpmorbit/Makefile \
                bpmprocess/Makefile \
                bpmrf/Makefile \
                bpmsimulation/Makefile \
                bpmdsp/Makefile \
                bpmwf/Makefile ])
AC_OUTPUT
