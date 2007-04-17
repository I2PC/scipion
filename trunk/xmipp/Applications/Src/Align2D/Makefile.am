bin_PROGRAMS = align2d
align2d_SOURCES = align2d.cc
align2d_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

align2d_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
