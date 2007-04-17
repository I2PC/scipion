bin_PROGRAMS = fourier_filter
fourier_filter_SOURCES = main.cpp
fourier_filter_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

fourier_filter_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
