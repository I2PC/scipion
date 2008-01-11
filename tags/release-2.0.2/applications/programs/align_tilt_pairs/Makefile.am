bin_PROGRAMS = align_tilt_pairs
align_tilt_pairs_SOURCES = main.cpp
align_tilt_pairs_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

align_tilt_pairs_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
