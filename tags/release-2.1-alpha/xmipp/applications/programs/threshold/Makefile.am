bin_PROGRAMS = threshold
threshold_SOURCES = main.cpp
threshold_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

threshold_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
