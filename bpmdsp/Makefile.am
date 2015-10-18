lib_LTLIBRARIES = libbpmdsp.la
libbpmdsp_la_SOURCES = create_filter.c \
                       apply_filter.c \
                       delete_filter.c \
                       print_filter.c \
                       calculate_filter_coefficients.c \
                       gaussian_filter_coeffs.c \
                       create_splane_representation.c \
                       create_resonator_representation.c \
                       normalise_filter.c \
                       print_filter_representation.c \
                       zplane_transform.c \
                       filter_step_response.c \
                       filter_impulse_response.c \
                       ddc.c fftsg.c norm_phase.c \
                       discrete_fourier_transforms.c


pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_dsp.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmdsp_la_LDFLAGS = -release $(BPM_LIB_VERSION)
