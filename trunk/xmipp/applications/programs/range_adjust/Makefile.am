bin_PROGRAMS = range_adjust
range_adjust_SOURCES = main.cpp
range_adjust_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

range_adjust_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
