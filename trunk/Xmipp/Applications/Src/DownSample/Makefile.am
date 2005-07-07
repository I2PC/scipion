bin_PROGRAMS = downsample
downsample_SOURCES = downsample.cc xvsmooth.c
downsample_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

downsample_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

EXTRA_DIST = xvsmooth.h
