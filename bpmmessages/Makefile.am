lib_LTLIBRARIES = libbpmmessages.la 
libbpmmessages_la_SOURCES = bpm_error.c \
                            bpm_warning.c

pkgincludedir=$(includedir)/bpm
pkginclude_HEADERS = bpm_messages.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libbpmmessages_la_LDFLAGS = -release $(BPM_LIB_VERSION)

