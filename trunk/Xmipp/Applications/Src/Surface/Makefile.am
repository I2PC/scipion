bin_PROGRAMS = surface
surface_SOURCES = surface.cc
surface_LDADD =$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

surface_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/Lib

