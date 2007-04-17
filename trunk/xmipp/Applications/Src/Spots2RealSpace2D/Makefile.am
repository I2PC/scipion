bin_PROGRAMS = crystal_aph2img
crystal_aph2img_SOURCES = spots2realspace.cc
crystal_aph2img_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
crystal_aph2img_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
