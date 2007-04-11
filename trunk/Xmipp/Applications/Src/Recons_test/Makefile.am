bin_PROGRAMS = phantom_test_reconstruction
phantom_test_reconstruction_SOURCES = recons_test.cc
phantom_test_reconstruction_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

phantom_test_reconstruction_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
