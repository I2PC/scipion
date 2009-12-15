bin_PROGRAMS = ml_align2D_combine
ml_align2D_combine_SOURCES = main.cpp
ml_align2D_combine_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

ml_align2D_combine_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
