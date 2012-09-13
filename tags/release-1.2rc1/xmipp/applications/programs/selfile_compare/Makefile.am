bin_PROGRAMS = selfile_compare
selfile_compare_SOURCES = main.cpp
selfile_compare_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

selfile_compare_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
