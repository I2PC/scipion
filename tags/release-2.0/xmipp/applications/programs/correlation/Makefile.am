bin_PROGRAMS = correlation
correlation_SOURCES = main.cpp
correlation_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

correlation_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
