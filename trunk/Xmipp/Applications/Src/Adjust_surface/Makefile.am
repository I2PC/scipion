bin_PROGRAMS = adjust_surface
adjust_surface_SOURCES = adjust_surface.cc
adjust_surface_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippReconsINRIA_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

adjust_surface_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
