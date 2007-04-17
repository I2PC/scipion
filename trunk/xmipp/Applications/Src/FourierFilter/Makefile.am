bin_PROGRAMS = fourier_filter
fourier_filter_SOURCES = fourierfilter.cc
fourier_filter_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
fourier_filter_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
