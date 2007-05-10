bin_PROGRAMS = normalize
normalize_SOURCES = main.cpp
normalize_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

normalize_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries

