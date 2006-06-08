bin_PROGRAMS = EM2Angles
EM2Angles_SOURCES = EM2Angles.cc
EM2Angles_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

EM2Angles_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
