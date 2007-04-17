bin_PROGRAMS = crystal_skew
crystal_skew_SOURCES = main.cpp
crystal_skew_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

crystal_skew_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
