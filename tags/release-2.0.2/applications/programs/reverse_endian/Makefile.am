bin_PROGRAMS = reverse_endian
reverse_endian_SOURCES = main.cpp
reverse_endian_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

reverse_endian_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
