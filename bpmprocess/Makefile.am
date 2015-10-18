lib_LTLIBRARIES = libbpmprocess.la 
libbpmprocess_la_SOURCES = process_diode.c \
                           process_monopole.c \
                           process_waveform.c \
                           postprocess_waveform.c \
                           process_caltone.c \
                           process_dipole.c \
                           ddc_waveform.c \
                           ddc_sample_waveform.c \
                           fft_waveform.c \
                           fit_diodepulse.c \
                           fit_fft.c \
                           fit_waveform.c \
                           get_IQ.c \
                           get_pedestal.c \
                           get_pos.c \
                           get_t0.c \
                           get_slope.c \
                           check_saturation.c \
                           correct_gain.c \
                           downmix_waveform.c


pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_process.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)


libbpmprocess_la_LDFLAGS = -release $(BPM_LIB_VERSION)
