bin_PROGRAMS = microscope
microscope_SOURCES = microscope.cc
microscope_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
INCLUDES = -I$(top_srcdir)/Lib

