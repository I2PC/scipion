bin_PROGRAMS = radon_transform
radon_transform_SOURCES = radon_transform.cc
radon_transform_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	

INCLUDES = -I$(top_srcdir)/Lib

