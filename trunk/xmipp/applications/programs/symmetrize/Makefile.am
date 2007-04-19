bin_PROGRAMS = symmetrize
symmetrize_SOURCES = main.cpp
symmetrize_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

symmetrize_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
