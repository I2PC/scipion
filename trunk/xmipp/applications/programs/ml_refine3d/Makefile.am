bin_PROGRAMS = ml_refine3d
ml_refine3d_SOURCES = main.cpp
ml_refine3d_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

ml_refine3d_LDFLAGS = @EXEC_TYPE@
INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
