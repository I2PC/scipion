bin_PROGRAMS = correlation
correlation_SOURCES = main.cpp
correlation_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la

correlation_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
