bin_PROGRAMS = IFFT
IFFT_SOURCES = IFFT.cc
IFFT_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la

IFFT_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
