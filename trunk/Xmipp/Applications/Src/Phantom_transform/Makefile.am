bin_PROGRAMS = phantom_transform
phantom_transform_SOURCES = phantom_transform.cc
phantom_transform_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

phantom_transform_LDFLAGS = @EXEC_TYPE@
	
INCLUDES = -I$(top_srcdir)/Lib


