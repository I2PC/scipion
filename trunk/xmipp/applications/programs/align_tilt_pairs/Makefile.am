bin_PROGRAMS = align_tilt_pairs
align_tilt_pairs_SOURCES = centilt.cc
align_tilt_pairs_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

align_tilt_pairs_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
