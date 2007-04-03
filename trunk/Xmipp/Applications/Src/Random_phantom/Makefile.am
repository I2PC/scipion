bin_PROGRAMS = create_random_phantom
create_random_phantom_SOURCES = random_phantom.cc
create_random_phantom_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
create_random_phantom_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
