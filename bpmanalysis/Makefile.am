lib_LTLIBRARIES = libbpmanalysis.la 
libbpmanalysis_la_SOURCES = ana_set_cutfn.c \
                            ana_def_cutfn.c \
                            ana_get_svd_coeffs.c \
                            ana_compute_residual.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_analysis.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmanalysis_la_LDFLAGS = -release $(BPM_LIB_VERSION)

