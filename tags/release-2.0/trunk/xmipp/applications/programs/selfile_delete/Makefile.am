bin_PROGRAMS = selfile_delete
selfile_delete_SOURCES = main.cpp
selfile_delete_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

selfile_delete_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
