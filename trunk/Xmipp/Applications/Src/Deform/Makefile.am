bin_PROGRAMS = deform
deform_SOURCES = deform.cc
deform_LDADD =$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
 	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
	
INCLUDES = -I$(top_srcdir)/Lib
