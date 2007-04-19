bin_PROGRAMS = ml_align2d
ml_align2d_SOURCES = main.cpp
ml_align2d_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

ml_align2d_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
