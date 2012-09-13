bin_PROGRAMS = histogram
histogram_SOURCES = main.cpp
histogram_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

histogram_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries


