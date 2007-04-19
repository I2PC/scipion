bin_PROGRAMS = adjust_volume_grey_levels
adjust_volume_grey_levels_SOURCES = main.cpp
adjust_volume_grey_levels_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

adjust_volume_grey_levels_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
