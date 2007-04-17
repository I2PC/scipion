bin_PROGRAMS = micrograph_downsample
micrograph_downsample_SOURCES = main.cpp xvsmooth.cpp
micrograph_downsample_LDADD = $(top_srcdir)/libraries/data/libXmippData.la
micrograph_downsample_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
EXTRA_DIST = xvsmooth.h
