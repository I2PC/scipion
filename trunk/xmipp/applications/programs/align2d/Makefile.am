bin_PROGRAMS = align2d
align2d_SOURCES = main.cpp
align2d_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

align2d_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
