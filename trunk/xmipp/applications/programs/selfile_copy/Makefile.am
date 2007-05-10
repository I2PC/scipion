bin_PROGRAMS = selfile_copy
selfile_copy_SOURCES = main.cpp
selfile_copy_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

selfile_copy_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
