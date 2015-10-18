lib_LTLIBRARIES = libbpmnr.la 
libbpmnr_la_SOURCES = nr_fit.c \
                      nr_four1.c \
                      nr_gammln.c \
                      nr_gammq.c \
                      nr_gcf.c \
                      nr_gser.c \
                      nr_realft.c \
                      nr_checks.c \
                      nr_rangauss.c \
                      nr_ran1.c \
                      nr_ranuniform.c \
                      nr_seed.c \
                      nr_levmar.c \
                      nr_median.c \
                      nr_select.c \
                      gsl_matrix.c \
                      gsl_vector.c \
                      gsl_linalg.c \
                      gsl_block.c \
                      gsl_blas.c \
                      nr_complex.c \
                      nr_sinc.c \
                      nr_quadinterpol.c \
                      dround.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_nr.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmnr_la_LDFLAGS = -release $(BPM_LIB_VERSION)

