bin_PROGRAMS = spots2realspace
spots2realspace_SOURCES = spots2realspace.cc
spots2realspace_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
	
spots2realspace_LDFLAGS = -all-static

INCLUDES = -I$(top_srcdir)/Lib
