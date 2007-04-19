bin_PROGRAMS = shift
shift_SOURCES = main.cpp
shift_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

shift_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
