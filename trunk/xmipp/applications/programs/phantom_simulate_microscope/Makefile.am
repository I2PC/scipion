bin_PROGRAMS = phantom_simulate_microscope
phantom_simulate_microscope_SOURCES = main.cpp
phantom_simulate_microscope_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

phantom_simulate_microscope_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries

