bin_PROGRAMS = crystal_euler2mrc
crystal_euler2mrc_SOURCES = EM2Angles.cc
crystal_euler2mrc_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

crystal_euler2mrc_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
