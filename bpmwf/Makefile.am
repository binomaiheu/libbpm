lib_LTLIBRARIES = libbpmwf.la
libbpmwf_la_SOURCES = doublewf.c \
                      intwf.c \
                      complexwf.c \
                      wfstats.c \
                      time_to_sample.c \
                      sample_to_time.c \
                      freq_to_sample.c \
                      sample_to_freq.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_wf.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmwf_la_LDFLAGS = -release $(BPM_LIB_VERSION)
