bin_PROGRAMS = sort_by_statistics
sort_by_statistics_SOURCES = main.cpp
sort_by_statistics_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

sort_by_statistics_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
