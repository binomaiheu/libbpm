lib_LTLIBRARIES = libbpmrf.la 
libbpmrf_la_SOURCES = rf_addLO.c \
                      rf_amplify.c \
		      rf_amplify_complex.c \
		      rf_phase_shifter.c \
                      rf_mixer.c \
                      rf_rectify.c \
		      rf_setup.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_rf.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)


libbpmrf_la_LDFLAGS = -release $(BPM_LIB_VERSION)
