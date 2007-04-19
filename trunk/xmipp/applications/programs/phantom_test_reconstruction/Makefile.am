bin_PROGRAMS = phantom_test_reconstruction
phantom_test_reconstruction_SOURCES = main.cpp
phantom_test_reconstruction_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons_Interface.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

phantom_test_reconstruction_LDFLAGS = @EXEC_TYPE@

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
