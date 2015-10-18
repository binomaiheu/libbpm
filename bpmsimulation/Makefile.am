lib_LTLIBRARIES = libbpmsimulation.la 
libbpmsimulation_la_SOURCES = digitise.c \
                              get_mode_response.c \
                              generate_bpmsignal.c \
                              add_mode_response.c \
                              get_mode_amplitude.c \
			      set_temp.c \
			      set_time.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_simulation.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)


libbpmsimulation_la_LDFLAGS = -release $(BPM_LIB_VERSION)
