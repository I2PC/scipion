bin_PROGRAMS = average
average_SOURCES = main.cpp
average_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

average_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
