bin_PROGRAMS = phantom_simulate_microscope
phantom_simulate_microscope_SOURCES = microscope.cc
phantom_simulate_microscope_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
phantom_simulate_microscope_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib

