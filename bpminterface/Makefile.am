lib_LTLIBRARIES = libbpminterface.la 
libbpminterface_la_SOURCES = bpm_verbose.c bpm_evtnum.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_interface.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)


libbpminterface_la_LDFLAGS = -release $(BPM_LIB_VERSION)
