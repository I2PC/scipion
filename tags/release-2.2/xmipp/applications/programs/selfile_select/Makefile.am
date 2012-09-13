bin_PROGRAMS = selfile_select
selfile_select_SOURCES = main.cpp
selfile_select_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

selfile_select_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
