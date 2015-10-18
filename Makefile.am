# automake 
# makefile template
#
# Bino Maiheu, 2006-2007, University College London

include $(top_srcdir)/aminclude.am

SUBDIRS = bpm \
          bpmanalysis \
          bpmcalibration \
          bpminterface \
          bpmmessages \
          bpmnr \
          bpmorbit \
          bpmprocess \
          bpmrf \
          bpmsimulation \
          bpmdsp \
          bpmwf

SUBLIBS = bpmanalysis/libbpmanalysis.la \
          bpmcalibration/libbpmcalibration.la \
          bpminterface/libbpminterface.la \
          bpmmessages/libbpmmessages.la \
          bpmnr/libbpmnr.la \
          bpmorbit/libbpmorbit.la \
          bpmprocess/libbpmprocess.la \
          bpmrf/libbpmrf.la \
          bpmsimulation/libbpmsimulation.la \
          bpmdsp/libbpmdsp.la \
          bpmwf/libbpmwf.la

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_defs.h bpm_version.h bpm_units.h

bin_SCRIPTS = libbpm-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= libbpm.pc

EXTRA_DIST = autogen.sh libbpm-config.in libbpm.pc.in configure.ac \
             acinclude.m4 aminclude.am doxygen.cfg \
             AUTHORS ChangeLog NEWS README INSTALL COPYING \
             documentation \
             bpm_version.h.in \
             examples/README examples/bpmtest.c examples/online_proc.c

lib_LTLIBRARIES   = libbpm.la
libbpm_la_SOURCES = version.c
libbpm_la_LIBADD  = $(SUBLIBS)
libbpm_la_LDFLAGS = -release $(BPM_LIB_VERSION)

bin_PROGRAMS = bpmtest
bpmtest_SOURCES = examples/bpmtest.c
bpmtest_LDADD = $(SUBLIBS) -lm 



# m4datadir = $(datadir)/aclocal
# m4data_DATA = bpm.m4

