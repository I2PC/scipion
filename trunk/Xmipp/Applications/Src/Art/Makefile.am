bin_PROGRAMS = Art
Art_SOURCES = art.cc
Art_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterfaceVTK_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
	
Art_LDFLAGS = -all-static

INCLUDES = -I$(top_srcdir)/Lib
