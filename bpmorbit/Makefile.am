lib_LTLIBRARIES = libbpmorbit.la 
libbpmorbit_la_SOURCES = get_bpmhit.c \
                         get_bend.c \
                         vm.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_orbit.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmorbit_la_LDFLAGS = -release $(BPM_LIB_VERSION)
