bin_PROGRAMS = fourierfilter
fourierfilter_SOURCES = fourierfilter.cc
fourierfilter_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
fourierfilter_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

