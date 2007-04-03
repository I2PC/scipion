bin_PROGRAMS = phantom_create
phantom_create_SOURCES = create_phantom.cc
phantom_create_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
phantom_create_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
