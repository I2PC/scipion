bin_PROGRAMS = random_phantom
random_phantom_SOURCES = random_phantom.cc
random_phantom_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
INCLUDES = -I$(top_srcdir)/Lib
