bin_PROGRAMS = normalize
normalize_SOURCES = normalize.cc
normalize_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

normalize_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/Lib

