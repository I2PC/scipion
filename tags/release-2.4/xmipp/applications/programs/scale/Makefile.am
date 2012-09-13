bin_PROGRAMS = scale
scale_SOURCES = main.cpp
scale_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

scale_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
