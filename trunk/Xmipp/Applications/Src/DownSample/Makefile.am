bin_PROGRAMS = micrograph_downsample
micrograph_downsample_SOURCES = downsample.cc xvsmooth.cc
micrograph_downsample_LDADD = \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
micrograph_downsample_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
EXTRA_DIST = xvsmooth.h
