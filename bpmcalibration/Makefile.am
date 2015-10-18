lib_LTLIBRARIES = libbpmcalibration.la 
libbpmcalibration_la_SOURCES = calibrate.c \
                               setup_calibration.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_calibration.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmcalibration_la_LDFLAGS = -release $(BPM_LIB_VERSION)
