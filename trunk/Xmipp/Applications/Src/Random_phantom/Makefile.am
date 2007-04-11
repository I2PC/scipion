bin_PROGRAMS = phantom_create_random
phantom_create_random_SOURCES = random_phantom.cc
phantom_create_random_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

phantom_create_randomm_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
