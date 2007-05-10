bin_PROGRAMS = header_print
header_print_SOURCES = main.cpp
header_print_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

header_print_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
