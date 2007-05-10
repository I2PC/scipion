bin_PROGRAMS = selfile_statistics
selfile_statistics_SOURCES = main.cpp
selfile_statistics_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

selfile_statistics_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
