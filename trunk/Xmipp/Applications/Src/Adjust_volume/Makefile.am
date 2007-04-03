bin_PROGRAMS = adjust_volume_grey_levels
adjust_volume_grey_levels_SOURCES = adjust_volume.cc
adjust_volume_grey_levels_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
adjust_volume_grey_levels_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
