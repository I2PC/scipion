bin_PROGRAMS = header_reset
header_reset_SOURCES = main.cpp
header_reset_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

header_reset_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
