bin_PROGRAMS = header_extract
header_extract_SOURCES = main.cpp
header_extract_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

header_extract_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
