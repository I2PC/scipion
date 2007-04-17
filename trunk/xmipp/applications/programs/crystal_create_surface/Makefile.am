bin_PROGRAMS = create_surface
create_surface_SOURCES = create_surface.cc
create_surface_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
                       $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
create_surface_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib


