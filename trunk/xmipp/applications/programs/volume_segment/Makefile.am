bin_PROGRAMS = volume_segment
volume_segment_SOURCES = main.cpp
volume_segment_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

volume_segment_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
